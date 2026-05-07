#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <unordered_map>
#include <omp.h>

#include "utils.h"
#include "known_ranks.h"
#include "parameters/flip_parameters.h"
#include "parameters/meta_parameters.h"

template <typename Scheme>
class MetaFlipGraph {
    size_t count;
    std::string outputPath;
    int threads;
    FlipParameters flipParameters;
    MetaParameters metaParameters;
    int seed;
    size_t topCount;
    std::string format;

    std::vector<Scheme> schemes;
    std::vector<Scheme> schemesBest;
    std::vector<size_t> flips;
    std::vector<size_t> iterations;
    std::vector<size_t> plusIterations;
    std::vector<int> bestRanks;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
    std::uniform_int_distribution<size_t> plusDistribution;

    std::unordered_map<std::string, std::vector<Scheme>> dimension2improvements;
    std::unordered_map<std::string, int> dimension2bestRank;
    std::unordered_map<std::string, int> dimension2knownRank;
    std::unordered_map<std::string, std::vector<int>> dimension2indices;
    std::vector<std::string> dimensionsInitial;
    std::vector<std::string> dimensions;
public:
    MetaFlipGraph(size_t count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const MetaParameters &metaParameters, int seed, size_t topCount, const std::string &format);

    bool initializeNaive(int n1, int n2, int n3);
    bool initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness);

    void initializeKnownRanks(const std::string &ring);
    void run();
private:
    void initialize();
    void flipIteration();
    void metaIteration();

    void updateBest(size_t iteration);
    void updateRanks(int iteration, bool save);
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator);
    void meta(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::mt19937 &generator);

    bool metaDefault(Scheme &scheme, std::mt19937 &generator);
    bool metaProjections(Scheme &scheme, std::mt19937 &generator);
    bool metaExtensions(Scheme &scheme, std::mt19937 &generator);

    void updateIndices();
    bool compare(int index1, int index2) const;
    std::string getSavePath(const Scheme &scheme, int iteration, const std::string path) const;
    std::string sortedDimension(const Scheme &scheme) const;

    void saveScheme(const Scheme &scheme, const std::string &path) const;
    bool resetToNormalScheme(Scheme &scheme, std::mt19937 &generator);
};

template <typename Scheme>
MetaFlipGraph<Scheme>::MetaFlipGraph(size_t count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const MetaParameters &metaParameters, int seed, size_t topCount, const std::string &format) : uniform(0.0, 1.0), plusDistribution(flipParameters.minPlusIterations, flipParameters.maxPlusIterations) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, (int) count);

    this->flipParameters = flipParameters;
    this->metaParameters = metaParameters;

    this->seed = seed;
    this->topCount = std::min(topCount, count);
    this->format = format;

    generators = initRandomGenerators(seed, threads);

    schemes.resize(count);
    schemesBest.resize(count);
    bestRanks.resize(count);
    flips.resize(count);
    iterations.resize(count);
    plusIterations.resize(count);
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::initializeNaive(int n1, int n2, int n3) {
    std::cout << "Start initializing with naive " << n1 << "x" << n2 << "x" << n3 << " schemes" << std::endl;

    if (!schemes[0].initializeNaive(n1, n2, n3))
        return false;

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 1; i < count; i++)
        schemes[i].initializeNaive(n1, n2, n3);

    std::string dimension = sortedDimension(schemes[0]);
    dimension2improvements.clear();
    dimension2improvements[dimension].push_back(Scheme(schemes[0]));
    dimensionsInitial.push_back(dimension);
    return true;
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness) {
    std::ifstream f(path);
    if (!f) {
        std::cout << "error: unable to open file \"" << path << "\"" << std::endl;
        return false;
    }

    bool valid = true;
    size_t schemesCount = 1;

    if (multiple)
        f >> schemesCount;

    std::cout << "Start reading " << std::min(count, schemesCount) << " / " << schemesCount << " schemes from \"" << path << "\"" << std::endl;

    for (size_t i = 0; i < count && i < schemesCount && valid; i++) {
        valid = schemes[i].read(f, checkCorrectness);

        if (!valid)
            std::cout << "error: invalid scheme " << (i + 1) << " in the file \"" << path << "\"" << std::endl;
    }

    f.close();

    if (!valid)
        return false;

    dimension2improvements.clear();
    dimensionsInitial.clear();

    for (size_t i = 0; i < count && i < schemesCount; i++) {
        std::string dimension = sortedDimension(schemes[i]);
        dimension2improvements[dimension].push_back(Scheme(schemes[i]));

        if (dimension2knownRank.find(dimension) == dimension2knownRank.end() || schemes[i].getRank() < dimension2knownRank.at(dimension))
            dimension2knownRank[dimension] = schemes[i].getRank();
    }

    for (auto &pair : dimension2improvements) {
        std::sort(pair.second.begin(), pair.second.end(), [](const Scheme& s1, const Scheme &s2) { return s1.getRank() > s2.getRank(); });
        dimensionsInitial.push_back(pair.first);
    }

    #pragma omp parallel for num_threads(threads)
    for (size_t i = schemesCount; i < count; i++)
        schemes[i].copy(schemes[i % schemesCount]);

    return true;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initializeKnownRanks(const std::string &ring) {
    dimension2knownRank = KNOWN_RANKS.at(ring);
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::run() {
    initialize();

    auto startTime = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    for (size_t iteration = 0; 1; iteration++) {
        flipIteration();
        updateBest(iteration);
        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);

        report(iteration + 1, startTime, elapsedTimes);

        t1 = std::chrono::high_resolution_clock::now();

        if (metaParameters.probability > 0) {
            metaIteration();
            updateRanks(iteration, true);
        }
    }
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initialize() {
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < count; i++) {
        int rank = schemes[i].getRank();

        schemesBest[i].copy(schemes[i]);
        bestRanks[i] = rank;
        flips[i] = 0;
        iterations[i] = 0;
        plusIterations[i] = plusDistribution(generators[omp_get_thread_num()]);

        dimension2bestRank[sortedDimension(schemes[i])] = rank;
    }

    updateRanks(0, false);
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::flipIteration()  {
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < count; i++)
        randomWalk(schemes[i], schemesBest[i], flips[i], iterations[i], plusIterations[i], bestRanks[i], generators[omp_get_thread_num()]);
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::metaIteration() {
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < count; i++)
        meta(schemes[i], flips[i], iterations[i], plusIterations[i], generators[omp_get_thread_num()]);
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::updateBest(size_t iteration) {
    updateIndices();

    for (const auto &pair : dimension2indices) {
        int top = pair.second[0];
        int bestRank = dimension2bestRank[pair.first];

        if (bestRanks[top] >= bestRank)
            continue;

        auto known = dimension2knownRank.find(pair.first);
        if (known != dimension2knownRank.end() && bestRanks[top] >= known->second)
            continue;

        if (!schemesBest[top].validateParallel()) {
            std::cout << "error: unable to save scheme " << schemesBest[top].getDimension() << " - it is invalid" << std::endl;
            exit(-1);
        }

        std::string path = getSavePath(schemesBest[top], iteration, outputPath);
        saveScheme(schemesBest[top], path);
        dimension2improvements[pair.first].push_back(Scheme(schemesBest[top]));

        std::cout << "Rank of " << pair.first << " was improved from " << bestRank << " to " << bestRanks[top] << ", scheme was saved to \"" << path << "." << format << "\"" << std::endl;
        dimension2bestRank[pair.first] = bestRanks[top];
    }
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::updateRanks(int iteration, bool save) {
    std::vector<std::string> newDimensions(count);
    std::unordered_map<std::string, int> dimension2bestIndex;

    for (size_t i = 0; i < count; i++) {
        std::string dimension = sortedDimension(schemes[i]);
        int rank = schemes[i].getRank();
        newDimensions[i] = dimension;

        auto result = dimension2bestRank.find(dimension);
        if (result == dimension2bestRank.end() || rank < result->second) {
            dimension2bestRank[dimension] = rank;
            dimension2bestIndex[dimension] = i;
        }

        if (sortedDimension(schemesBest[i]) != dimension)
            schemesBest[i].copy(schemes[i]);
    }

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < count; i++)
        bestRanks[i] = dimension2bestRank[newDimensions[i]];

    for (const auto &pair : dimension2bestIndex) {
        if (dimension2improvements.find(pair.first) == dimension2improvements.end())
            dimension2improvements[pair.first].push_back(Scheme(schemes[pair.second]));

        if (!save)
            continue;

        auto known = dimension2knownRank.find(pair.first);
        if (known != dimension2knownRank.end() && bestRanks[pair.second] >= known->second)
            continue;

        if (!schemes[pair.second].validateParallel()) {
            std::cout << "error: unable to save scheme " << schemes[pair.second].getDimension() << " - it is invalid" << std::endl;
            exit(-1);
        }

        std::string path = getSavePath(schemes[pair.second], iteration, outputPath);
        saveScheme(schemes[pair.second], path);
        std::cout << "Rank of " << pair.first << " was improved to " << bestRanks[pair.second] << ", scheme was saved to \"" << path << "." << format << "\"" << std::endl;
    }
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    std::cout << "+---------------------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| " << std::left;
    std::cout << "ring: " << std::setw(21) << schemes[0].getRing() << "   ";
    std::cout << "count: " << std::setw(20) << count << "   ";
    std::cout << std::right << std::setw(43) << ("iteration: " + std::to_string(iteration));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "seed: " << std::setw(21) << seed << "   ";
    std::cout << "threads: " << std::setw(18) << threads << "   ";
    std::cout << std::right << std::setw(43) << ("elapsed: " + prettyTime(elapsed));
    std::cout << " |" << std::endl;

    bool improved = false;
    for (const auto &pair : dimension2bestRank) {
        std::string dimension = pair.first;
        int rank = pair.second;

        auto known = dimension2knownRank.find(dimension);
        if (known == dimension2knownRank.end() || rank >= known->second)
            continue;

        if (!improved) {
            std::cout << "+=========================================================================================================+" << std::endl;
            std::cout << "| Improvements:                                                                                           |" << std::endl;
            std::cout << "| +-----------+------------+---------------+                                                              |" << std::endl;
            std::cout << "| | dimension | known rank | improved rank |                                                              |" << std::endl;
            std::cout << "| +-----------+------------+---------------+                                                              |" << std::endl;
            improved = true;
        }

        std::cout << "| | " << std::setw(9) << dimension << " | " << std::setw(10) << known->second << " | " << std::setw(13) << rank << " |                                                              |" << std::endl;
    }

    if (improved) {
        std::cout << "| +-----------+------------+---------------+                                                              |" << std::endl;
        std::cout << "|                                                                                                         |" << std::endl;
    }

    std::cout << "+=========================================================================================================+" << std::endl;
    std::cout << "| runner |     scheme size     | scheme rank |   naive    |            |        flips        |    plus    |" << std::endl;
    std::cout << "|   id   |  sorted  |   real   | best | curr | complexity | iterations |  count  | available | iterations |" << std::endl;
    std::cout << "+--------+----------+----------+------+------+------------+------------+---------+-----------+------------+" << std::endl;
    std::cout << std::right;

    for (const auto &dimension : dimensions) {
        const std::vector<int> &indices = dimension2indices.at(dimension);

        for (size_t i = 0; i < topCount && i < indices.size(); i++) {
            int runner = indices[i];

            std::cout << "| ";
            std::cout << std::setw(6) << runner << " | ";
            std::cout << std::setw(8) << dimension << " | ";
            std::cout << std::setw(8) << schemes[runner].getDimension() << " | ";
            std::cout << std::setw(4) << bestRanks[runner] << " | ";
            std::cout << std::setw(4) << schemes[runner].getRank() << " | ";
            std::cout << std::setw(10) << schemes[runner].getComplexity() << " | ";
            std::cout << std::setw(10) << prettyInt(iterations[runner]) << " | ";
            std::cout << std::setw(7) << prettyInt(flips[runner]) << " | ";
            std::cout << std::setw(9) << schemes[runner].getAvailableFlips() << " | ";
            std::cout << std::setw(10) << prettyInt(plusIterations[runner]) << " |";
            std::cout << std::endl;
        }

        std::cout << "+--------+----------+----------+------+------+------------+------------+---------+-----------+------------+" << std::endl;
    }

    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator) {
    plusIterations = plusDistribution(generator);

    for (size_t iteration = 0; iteration < flipParameters.flipIterations; iteration++) {
        int prevRank = scheme.getRank();

        if (!scheme.tryFlip(generator)) {
            if (scheme.tryExpand(generator))
                flipsCount = 0;

            continue;
        }

        if (flipParameters.reduceProbability && uniform(generator) < flipParameters.reduceProbability && scheme.tryReduce())
            flipsCount = 0;

        if (flipParameters.sandwichingProbability && uniform(generator) < flipParameters.sandwichingProbability)
            scheme.trySandwiching(generator);

        int rank = scheme.getRank();
        if (rank < prevRank)
            flipsCount = 0;

        flipsCount++;
        iterationsCount++;

        if (rank < bestRank) {
            schemeBest.copy(scheme);
            bestRank = rank;
            iterationsCount = 0;
        }

        if (flipsCount >= plusIterations && rank < bestRank + flipParameters.plusDiff && scheme.tryExpand(generator))
            flipsCount = 0;

        if (iterationsCount >= flipParameters.resetIterations) {
            std::string dimension = sortedDimension(scheme);
            Scheme &initial = dimension2improvements[dimension][generator() % dimension2improvements[dimension].size()];

            scheme.copy(initial);
            schemeBest.copy(initial);
            bestRank = initial.getRank();
            flipsCount = 0;
            iterationsCount = 0;
            plusIterations = plusDistribution(generator);
        }
    }
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::meta(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::mt19937 &generator) {
    if (uniform(generator) > metaParameters.probability)
        return;

    bool resized = false;

    if (metaParameters.strategy == "default") {
        resized = metaDefault(scheme, generator);
    }
    else if (metaParameters.strategy == "proj") {
        resized = metaProjections(scheme, generator);
    }
    else if (metaParameters.strategy == "ext") {
        resized = metaExtensions(scheme, generator);
    }

    if (!resized)
        return;

    flipsCount = 0;
    iterationsCount = 0;
    plusIterations = plusDistribution(generator);
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::metaDefault(Scheme &scheme, std::mt19937 &generator) {
    bool reset = resetToNormalScheme(scheme, generator);

    if (uniform(generator) < 0.5)
        scheme.swapSizes(generator);

    int index = generator() % schemesBest.size();

    if (scheme.tryMerge(schemesBest[index], generator, metaParameters.maxDimension, metaParameters.maxRank))
        return true;

    if (uniform(generator) < 0.5)
        return scheme.tryProject(generator, metaParameters.minDimension) || reset;

    return scheme.tryExtend(generator, metaParameters.maxDimension, metaParameters.maxRank) || reset;
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::metaProjections(Scheme &scheme, std::mt19937 &generator) {
    bool reset = resetToNormalScheme(scheme, generator);
    return scheme.tryProject(generator, metaParameters.minDimension) || reset;
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::metaExtensions(Scheme &scheme, std::mt19937 &generator) {
    bool reset = resetToNormalScheme(scheme, generator);
    return scheme.tryExtend(generator, metaParameters.maxDimension, metaParameters.maxRank) || reset;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::updateIndices() {
    dimension2indices.clear();
    dimensions.clear();

    for (size_t i = 0; i < count; i++)
        dimension2indices[sortedDimension(schemes[i])].push_back(i);

    for (auto &pair : dimension2indices) {
        size_t partialCount = std::min(topCount, pair.second.size());
        std::partial_sort(pair.second.begin(), pair.second.begin() + partialCount, pair.second.end(), [this](int index1, int index2) {
            return compare(index1, index2);
        });
    }

    for (const auto &pair : dimension2indices)
        dimensions.push_back(pair.first);

    std::sort(dimensions.begin(), dimensions.end(), [this](std::string &dimension1, std::string &dimension2) {
        const Scheme &scheme1 = schemes[dimension2indices.at(dimension1)[0]];
        const Scheme &scheme2 = schemes[dimension2indices.at(dimension2)[0]];
        int dimensions1[3];
        int dimensions2[3];

        for (int i = 0; i < 3; i++) {
            dimensions1[i] = scheme1.getDimension(i);
            dimensions2[i] = scheme2.getDimension(i);
        }

        std::sort(dimensions1, dimensions1 + 3);
        std::sort(dimensions2, dimensions2 + 3);

        for (int i = 0; i < 3; i++)
            if (dimensions1[i] != dimensions2[i])
                return dimensions1[i] < dimensions2[i];

        return false;
    });
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::compare(int index1, int index2) const {
    int bestRank1 = bestRanks[index1];
    int bestRank2 = bestRanks[index2];

    if (bestRank1 != bestRank2)
        return bestRank1 < bestRank2;

    int rank1 = schemes[index1].getRank();
    int rank2 = schemes[index2].getRank();

    if (rank1 != rank2)
        return rank1 < rank2;

    int complexity1 = schemes[index1].getComplexity();
    int complexity2 = schemes[index2].getComplexity();

    if (complexity1 != complexity2)
        return complexity1 < complexity2;

    return index1 < index2;
}

template <typename Scheme>
std::string MetaFlipGraph<Scheme>::getSavePath(const Scheme &scheme, int iteration, const std::string path) const {
    std::stringstream ss;
    ss << path << "/";
    ss << sortedDimension(scheme);
    ss << "_m" << scheme.getRank();
    ss << "_c" << scheme.getComplexity();
    ss << "_iteration" << iteration;
    ss << "_" << scheme.getDimension();
    ss << "_" << scheme.getRing();
    return ss.str();
}

template <typename Scheme>
std::string MetaFlipGraph<Scheme>::sortedDimension(const Scheme &scheme) const {
    int dimension[3] = {scheme.getDimension(0), scheme.getDimension(1), scheme.getDimension(2)};
    std::sort(dimension, dimension + 3);

    std::stringstream ss;
    ss << dimension[0] << "x" << dimension[1] << "x" << dimension[2];
    return ss.str();
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::saveScheme(const Scheme &scheme, const std::string &path) const {
    if (format == "json") {
        scheme.saveJson(path + ".json");
    }
    else if (format == "txt") {
        scheme.saveTxt(path + ".txt");
    }
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::resetToNormalScheme(Scheme &scheme, std::mt19937 &generator) {
    auto rank = dimension2knownRank.find(sortedDimension(scheme));

    if (rank == dimension2knownRank.end() || scheme.getRank() > rank->second + metaParameters.maxRankDiff) {
        std::string dimension = dimensionsInitial[generator() % dimensionsInitial.size()];
        Scheme &initial = dimension2improvements[dimension][generator() % dimension2improvements[dimension].size()];
        scheme.copy(initial);
        return true;
    }

    return false;
}
