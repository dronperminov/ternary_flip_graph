#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <omp.h>

#include "utils.h"

template <typename Scheme>
class FlipGraph {
    int count;
    std::string outputPath;
    int threads;
    size_t flipIterations;
    size_t resetIterations;
    int plusDiff;
    double reduceProbability;
    double copyBestProbability;
    int seed;
    int topCount;
    size_t maxImprovements;
    size_t improvementsIndex;
    std::string format;

    std::vector<Scheme> schemes;
    std::vector<Scheme> schemesBest;
    std::vector<Scheme> improvements;
    std::vector<size_t> flips;
    std::vector<size_t> iterations;
    std::vector<size_t> plusIterations;
    std::vector<int> bestRanks;
    std::vector<int> indices;
    int bestRank;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
    std::uniform_int_distribution<size_t> plusDistribution;
public:
    FlipGraph(int count, const std::string outputPath, int threads, size_t flipIterations, size_t minPlusIterations, size_t maxPlusIterations, size_t resetIterations, int plusDiff, double reduceProbability, double copyBestProbability, int seed, int topCount, size_t maxImprovements, const std::string &format);

    bool initializeNaive(int n1, int n2, int n3);
    bool initializeFromFile(const std::string &path, bool multiple);

    void run(int targetRank);
private:
    void resetImprovements();
    void addImprovement(const Scheme &scheme);

    void initialize();
    void runIteration();
    void updateBest(size_t iteration);
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator);

    bool compare(int index1, int index2) const;
    std::string getSavePath(const Scheme &scheme, int iteration, const std::string path) const;

    void saveScheme(const Scheme &scheme, const std::string &path) const;
};

template <typename Scheme>
FlipGraph<Scheme>::FlipGraph(int count, const std::string outputPath, int threads, size_t flipIterations, size_t minPlusIterations, size_t maxPlusIterations, size_t resetIterations, int plusDiff, double reduceProbability, double copyBestProbability, int seed, int topCount, size_t maxImprovements, const std::string &format) : uniform(0.0, 1.0), plusDistribution(minPlusIterations, maxPlusIterations) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, count);
    this->flipIterations = flipIterations;
    this->plusDiff = plusDiff;
    this->resetIterations = resetIterations;
    this->reduceProbability = reduceProbability;
    this->copyBestProbability = copyBestProbability;
    this->seed = seed;
    this->topCount = std::min(topCount, count);
    this->maxImprovements = maxImprovements;
    this->format = format;

    resetImprovements();

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

    schemes.resize(count);
    schemesBest.resize(count);
    bestRanks.resize(count);
    flips.resize(count);
    iterations.resize(count);
    plusIterations.resize(count);
    indices.resize(count);
}

template <typename Scheme>
bool FlipGraph<Scheme>::initializeNaive(int n1, int n2, int n3) {
    std::cout << "Start initializing with naive " << n1 << "x" << n2 << "x" << n3 << " schemes" << std::endl;

    if (!schemes[0].initializeNaive(n1, n2, n3))
        return false;

    #pragma omp parallel for num_threads(threads)
    for (int i = 1; i < count; i++)
        schemes[i].initializeNaive(n1, n2, n3);

    resetImprovements();
    addImprovement(schemes[0]);
    return true;
}

template <typename Scheme>
bool FlipGraph<Scheme>::initializeFromFile(const std::string &path, bool multiple) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    int schemesCount = 1;
    if (multiple)
        f >> schemesCount;

    std::cout << "Start reading " << std::min(schemesCount, count) << " / " << schemesCount << " schemes from \"" << path << "\"" << std::endl;

    bool valid = true;
    for (int i = 0; i < schemesCount && i < count && valid; i++)
        valid = schemes[i].read(f);

    f.close();

    if (!valid)
        return false;

    resetImprovements();
    for (int i = 0; i <schemesCount && i < count && improvements.size() < maxImprovements; i++)
        addImprovement(schemes[i]);

    #pragma omp parallel for num_threads(threads)
    for (int i = schemesCount; i < count; i++)
        schemes[i].copy(schemes[i % schemesCount]);

    return true;
}

template <typename Scheme>
void FlipGraph<Scheme>::run(int targetRank) {
    initialize();

    auto startTime = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    for (size_t iteration = 0; bestRank > targetRank; iteration++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        runIteration();
        updateBest(iteration);
        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);

        report(iteration + 1, startTime, elapsedTimes);
    }
}

template <typename Scheme>
void FlipGraph<Scheme>::resetImprovements() {
    improvements.clear();
    improvementsIndex = 0;
}

template <typename Scheme>
void FlipGraph<Scheme>::addImprovement(const Scheme &scheme) {
    if (improvements.size() < maxImprovements) {
        improvements.push_back(Scheme(scheme));
        improvementsIndex = 0;
    }
    else {
        improvements[improvementsIndex].copy(scheme);
        improvementsIndex = (improvementsIndex + 1) % maxImprovements;
    }
}

template <typename Scheme>
void FlipGraph<Scheme>::initialize() {
    bestRank = schemes[0].getRank();

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++) {
        int rank = schemes[i].getRank();
        if (rank < bestRank)
            bestRank = rank;

        schemesBest[i].copy(schemes[i]);
        bestRanks[i] = rank;
        flips[i] = 0;
        iterations[i] = 0;
        plusIterations[i] = plusDistribution(generators[omp_get_thread_num()]);
        indices[i] = i;
    }
}

template <typename Scheme>
void FlipGraph<Scheme>::runIteration()  {
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        randomWalk(schemes[i], schemesBest[i], flips[i], iterations[i], plusIterations[i], bestRanks[i], generators[omp_get_thread_num()]);
}

template <typename Scheme>
void FlipGraph<Scheme>::updateBest(size_t iteration) {
    std::partial_sort(indices.begin(), indices.begin() + topCount, indices.end(), [this](int index1, int index2) {
        return compare(index1, index2);
    });

    int top = indices[0];
    if (bestRanks[top] >= bestRank)
        return;

    if (!schemesBest[top].validate()) {
        std::cout << "Unable to save: scheme invalid" << std::endl;
        return;
    }

    std::string path = getSavePath(schemesBest[top], iteration, outputPath);
    saveScheme(schemesBest[top], path);
    addImprovement(schemesBest[top]);

    std::cout << "Rank was improved from " << bestRank << " to " << bestRanks[top] << ", scheme was saved to \"" << path << "." << format << "\"" << std::endl;
    bestRank = bestRanks[top];

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++) {
        if (uniform(generators[omp_get_thread_num()]) >= copyBestProbability || i == top)
            continue;

        schemes[i].copy(schemesBest[top]);
        schemesBest[i].copy(schemesBest[top]);
        iterations[i] = 0;
        flips[i] = 0;
    }
}

template <typename Scheme>
void FlipGraph<Scheme>::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    std::cout << std::left;
    std::cout << "+-----------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| ";
    std::cout << "dimension: " << std::setw(14) << schemes[indices[0]].getDimension() << "   ";
    std::cout << "seed: " << std::setw(20) << seed << "   ";
    std::cout << std::right << std::setw(24) << ("best rank: " + std::to_string(bestRank));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "threads: " << std::setw(16) << threads << "   ";
    std::cout << "flip iters: " << std::setw(14) << prettyInt(flipIterations) << "   ";
    std::cout << std::right << std::setw(24) << ("iteration: " + std::to_string(iteration));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "count: " << std::setw(18) << count << "   ";
    std::cout << "reset iters: " << std::setw(13) << prettyInt(resetIterations) << "   ";
    std::cout << std::right << std::setw(24) << ("elapsed: " + prettyTime(elapsed));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "ring: " << std::setw(19) << schemes[0].getRing() << "   ";
    std::cout << "plus diff: " << std::setw(15) << plusDiff << "   ";
    std::cout << std::right << std::setw(24) << ("improvements: " + std::to_string(improvements.size()) + " / " + std::to_string(maxImprovements));
    std::cout << " |" << std::endl;

    std::cout << "+===================================================================================+" << std::endl;
    std::cout << "| runner | scheme rank |   naive    |            |        flips        |    plus    |" << std::endl;
    std::cout << "|   id   | best | curr | complexity | iterations |  count  | available | iterations |" << std::endl;
    std::cout << "+--------+------+------+------------+------------+---------+-----------+------------+" << std::endl;
    std::cout << std::right;

    for (int i = 0; i < topCount; i++) {
        int runner = indices[i];
        std::cout << "| ";
        std::cout << std::setw(6) << runner << " | ";
        std::cout << std::setw(4) << schemesBest[runner].getRank() << " | ";
        std::cout << std::setw(4) << schemes[runner].getRank() << " | ";
        std::cout << std::setw(10) << schemes[runner].getComplexity() << " | ";
        std::cout << std::setw(10) << prettyInt(iterations[runner]) << " | ";
        std::cout << std::setw(7) << prettyInt(flips[runner]) << " | ";
        std::cout << std::setw(9) << schemes[runner].getAvailableFlips() << " | ";
        std::cout << std::setw(10) << prettyInt(plusIterations[runner]) << " |";
        std::cout << std::endl;
    }

    std::cout << "+--------+------+------+------------+------------+---------+-----------+------------+" << std::endl;
    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
void FlipGraph<Scheme>::randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator) {
    plusIterations = plusDistribution(generator);

    for (size_t iteration = 0; iteration < flipIterations; iteration++) {
        int prevRank = scheme.getRank();

        if (!scheme.tryFlip(generator)) {
            if (scheme.tryExpand(generator))
                flipsCount = 0;

            continue;
        }

        if (uniform(generator) < reduceProbability && scheme.tryReduce())
            flipsCount = 0;

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

        if (flipsCount >= plusIterations && rank < bestRank + plusDiff && scheme.tryExpand(generator))
            flipsCount = 0;

        if (iterationsCount >= resetIterations) {
            Scheme &initial = improvements[generator() % improvements.size()];
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
bool FlipGraph<Scheme>::compare(int index1, int index2) const {
    int bestRank1 = schemesBest[index1].getRank();
    int bestRank2 = schemesBest[index2].getRank();

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
std::string FlipGraph<Scheme>::getSavePath(const Scheme &scheme, int iteration, const std::string path) const {
    std::stringstream ss;
    ss << path << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_c" << scheme.getComplexity();
    ss << "_iteration" << iteration;
    ss << "_" << scheme.getRing();
    return ss.str();
}

template <typename Scheme>
void FlipGraph<Scheme>::saveScheme(const Scheme &scheme, const std::string &path) const {
    if (format == "json") {
        scheme.saveJson(path + ".json");
    }
    else if (format == "txt") {
        scheme.saveTxt(path + ".txt");
    }
}
