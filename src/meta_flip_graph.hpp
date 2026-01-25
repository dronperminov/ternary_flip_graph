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

template <typename Scheme>
class MetaFlipGraph {
    size_t count;
    std::string outputPath;
    int threads;
    size_t flipIterations;
    size_t resetIterations;
    int plusDiff;
    double sandwichingProbability;
    double reduceProbability;
    double resizeProbability;
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
    std::vector<std::string> dimensions;
public:
    MetaFlipGraph(size_t count, const std::string outputPath, int threads, size_t flipIterations, size_t minPlusIterations, size_t maxPlusIterations, size_t resetIterations, int plusDiff, double sandwichingProbability, double reduceProbability, double resizeProbability, int seed, size_t topCount, const std::string &format);

    bool initializeNaive(int n1, int n2, int n3);
    bool initializeFromFile(const std::string &path, bool multiple);
    void initializeBestTernaryRanks();
    void initializeBestBinaryRanks();

    void run();
private:
    void initialize();
    void runIteration();
    void resizeIteration();
    void updateBest(size_t iteration);
    void updateRanks(int iteration, bool save);
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator);
    void resize(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::mt19937 &generator);

    void updateIndices();
    bool compare(int index1, int index2) const;
    std::string getSavePath(const Scheme &scheme, int iteration, const std::string path) const;
    std::string sortedDimension(const Scheme &scheme) const;

    void saveScheme(const Scheme &scheme, const std::string &path) const;
};

template <typename Scheme>
MetaFlipGraph<Scheme>::MetaFlipGraph(size_t count, const std::string outputPath, int threads, size_t flipIterations, size_t minPlusIterations, size_t maxPlusIterations, size_t resetIterations, int plusDiff, double sandwichingProbability, double reduceProbability, double resizeProbability, int seed, size_t topCount, const std::string &format) : uniform(0.0, 1.0), plusDistribution(minPlusIterations, maxPlusIterations) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, (int) count);
    this->flipIterations = flipIterations;
    this->plusDiff = plusDiff;
    this->sandwichingProbability = sandwichingProbability;
    this->resetIterations = resetIterations;
    this->reduceProbability = reduceProbability;
    this->resizeProbability = resizeProbability;
    this->seed = seed;
    this->topCount = std::min(topCount, count);
    this->format = format;

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

    schemes.resize(count);
    schemesBest.resize(count);
    bestRanks.resize(count);
    flips.resize(count);
    iterations.resize(count);
    plusIterations.resize(count);

    dimension2knownRank = {
        {"2x2x2", 7}, {"2x2x3", 11}, {"2x2x4", 14}, {"2x2x5", 18}, {"2x2x6", 21}, {"2x2x7", 25}, {"2x2x8", 28}, {"2x2x9", 32}, {"2x2x10", 35}, {"2x2x11", 39}, {"2x2x12", 42}, {"2x2x13", 46}, {"2x2x14", 49}, {"2x2x15", 53}, {"2x2x16", 56},
        {"2x3x3", 15}, {"2x3x4", 20}, {"2x3x5", 25}, {"2x3x6", 30}, {"2x3x7", 35}, {"2x3x8", 40}, {"2x3x9", 45}, {"2x3x10", 50}, {"2x3x11", 55}, {"2x3x12", 60}, {"2x3x13", 65}, {"2x3x14", 70}, {"2x3x15", 75}, {"2x3x16", 80},
        {"2x4x4", 26}, {"2x4x5", 32}, {"2x4x6", 39}, {"2x4x7", 45}, {"2x4x8", 51}, {"2x4x9", 58}, {"2x4x10", 64}, {"2x4x11", 71}, {"2x4x12", 77}, {"2x4x13", 83}, {"2x4x14", 90}, {"2x4x15", 96}, {"2x4x16", 102},
        {"2x5x5", 40}, {"2x5x6", 47}, {"2x5x7", 55}, {"2x5x8", 63}, {"2x5x9", 72}, {"2x5x10", 79}, {"2x5x11", 87}, {"2x5x12", 94},
        {"2x6x6", 56}, {"2x6x7", 66}, {"2x6x8", 75}, {"2x6x9", 86}, {"2x6x10", 94},
        {"2x7x7", 76}, {"2x7x8", 88}, {"2x7x9", 99},
        {"2x8x8", 100},
        {"3x3x3", 23}, {"3x3x4", 29}, {"3x3x5", 36}, {"3x3x6", 40}, {"3x3x7", 49}, {"3x3x8", 55}, {"3x3x9", 63}, {"3x3x10", 69}, {"3x3x11", 76}, {"3x3x12", 80}, {"3x3x13", 89}, {"3x3x14", 95}, {"3x3x15", 103}, {"3x3x16", 109},
        {"3x4x4", 38}, {"3x4x5", 47}, {"3x4x6", 54}, {"3x4x7", 63}, {"3x4x8", 73}, {"3x4x9", 83}, {"3x4x10", 92}, {"3x4x11", 101}, {"3x4x12", 108}, {"3x4x13", 117}, {"3x4x14", 126}, {"3x4x15", 136}, {"3x4x16", 146},
        {"3x5x5", 58}, {"3x5x6", 68}, {"3x5x7", 79}, {"3x5x8", 90}, {"3x5x9", 104}, {"3x5x10", 115}, {"3x5x11", 126}, {"3x5x12", 136},
        {"3x6x6", 80}, {"3x6x7", 94}, {"3x6x8", 108}, {"3x6x9", 120}, {"3x6x10", 134},
        {"3x7x7", 111}, {"3x7x8", 126}, {"3x7x9", 142},
        {"3x8x8", 145},
        {"4x4x4", 48}, {"4x4x5", 61}, {"4x4x6", 73}, {"4x4x7", 85}, {"4x4x8", 96}, {"4x4x9", 104}, {"4x4x10", 120}, {"4x4x11", 130}, {"4x4x12", 142}, {"4x4x13", 152}, {"4x4x14", 165}, {"4x4x15", 177}, {"4x4x16", 189},
        {"4x5x5", 76}, {"4x5x6", 90}, {"4x5x7", 104}, {"4x5x8", 118}, {"4x5x9", 136}, {"4x5x10", 150}, {"4x5x11", 165}, {"4x5x12", 179},
        {"4x6x6", 105}, {"4x6x7", 123}, {"4x6x8", 140}, {"4x6x9", 159}, {"4x6x10", 175},
        {"4x7x7", 144}, {"4x7x8", 164}, {"4x7x9", 186},
        {"4x8x8", 182},
        {"5x5x5", 93}, {"5x5x6", 110}, {"5x5x7", 127}, {"5x5x8", 144}, {"5x5x9", 167}, {"5x5x10", 184}, {"5x5x11", 202}, {"5x5x12", 220},
        {"5x6x6", 130}, {"5x6x7", 150}, {"5x6x8", 170}, {"5x6x9", 197}, {"5x6x10", 217},
        {"5x7x7", 176}, {"5x7x8", 205}, {"5x7x9", 229},
        {"5x8x8", 230},
        {"6x6x6", 153}, {"6x6x7", 183}, {"6x6x8", 203}, {"6x6x9", 225}, {"6x6x10", 247},
        {"6x7x7", 215}, {"6x7x8", 239}, {"6x7x9", 268},
        {"6x8x8", 266},
        {"7x7x7", 249}, {"7x7x8", 277}, {"7x7x9", 315},
        {"7x8x8", 306},
        {"8x8x8", 336}
    };
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::initializeNaive(int n1, int n2, int n3) {
    std::cout << "Start initializing with naive " << n1 << "x" << n2 << "x" << n3 << " schemes" << std::endl;

    if (!schemes[0].initializeNaive(n1, n2, n3))
        return false;

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 1; i < count; i++)
        schemes[i].initializeNaive(n1, n2, n3);

    dimension2improvements.clear();
    dimension2improvements[sortedDimension(schemes[0])].push_back(Scheme(schemes[0]));
    return true;
}

template <typename Scheme>
bool MetaFlipGraph<Scheme>::initializeFromFile(const std::string &path, bool multiple) {
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
        valid = schemes[i].read(f);

        if (!valid)
            std::cout << "error: invalid scheme " << (i + 1) << " in the file \"" << path << "\"" << std::endl;
    }

    f.close();

    if (!valid)
        return false;

    dimension2improvements.clear();
    for (size_t i = 0; i < count && i < schemesCount; i++)
        dimension2improvements[sortedDimension(schemes[i])].push_back(Scheme(schemes[i]));

    #pragma omp parallel for num_threads(threads)
    for (size_t i = schemesCount; i < count; i++)
        schemes[i].copy(schemes[i % schemesCount]);

    return true;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initializeBestTernaryRanks() {
    dimension2knownRank["2x4x5"] = 33;
    dimension2knownRank["2x4x9"] = 59;
    dimension2knownRank["2x4x10"] = 65;
    dimension2knownRank["2x4x13"] = 84;
    dimension2knownRank["2x5x7"] = 57;
    dimension2knownRank["2x5x8"] = 65;
    dimension2knownRank["2x5x10"] = 80;
    dimension2knownRank["2x6x6"] = 57;
    dimension2knownRank["2x6x7"] = 68;
    dimension2knownRank["2x6x8"] = 77;
    dimension2knownRank["2x7x7"] = 77;
    dimension2knownRank["2x7x9"] = 102;
    dimension2knownRank["3x3x6"] = 42;
    dimension2knownRank["3x3x8"] = 57;
    dimension2knownRank["3x3x9"] = 64;
    dimension2knownRank["3x3x10"] = 71;
    dimension2knownRank["3x3x11"] = 78;
    dimension2knownRank["3x3x12"] = 84;
    dimension2knownRank["3x3x13"] = 91;
    dimension2knownRank["3x3x14"] = 98;
    dimension2knownRank["3x3x15"] = 106;
    dimension2knownRank["3x3x16"] = 113;
    dimension2knownRank["3x4x7"] = 64;
    dimension2knownRank["3x4x8"] = 74;
    dimension2knownRank["3x4x13"] = 118;
    dimension2knownRank["3x4x14"] = 128;
    dimension2knownRank["3x4x15"] = 137;
    dimension2knownRank["3x5x6"] = 70;
    dimension2knownRank["3x5x7"] = 83;
    dimension2knownRank["3x5x8"] = 94;
    dimension2knownRank["3x5x9"] = 105;
    dimension2knownRank["3x5x11"] = 128;
    dimension2knownRank["3x5x12"] = 140;
    dimension2knownRank["3x6x6"] = 83;
    dimension2knownRank["3x6x7"] = 96;
    dimension2knownRank["3x6x9"] = 124;
    dimension2knownRank["3x6x10"] = 137;
    dimension2knownRank["3x7x7"] = 113;
    dimension2knownRank["3x7x8"] = 128;
    dimension2knownRank["3x7x9"] = 145;
    dimension2knownRank["3x8x8"] = 148;
    dimension2knownRank["4x4x4"] = 49;
    dimension2knownRank["4x4x9"] = 110;
    dimension2knownRank["4x4x10"] = 122;
    dimension2knownRank["4x4x11"] = 134;
    dimension2knownRank["4x4x12"] = 145;
    dimension2knownRank["4x4x13"] = 157;
    dimension2knownRank["4x4x14"] = 169;
    dimension2knownRank["4x4x15"] = 181;
    dimension2knownRank["4x4x16"] = 192;
    dimension2knownRank["4x5x9"] = 137;
    dimension2knownRank["4x5x10"] = 151;
    dimension2knownRank["4x7x7"] = 145;
    dimension2knownRank["4x7x9"] = 187;
    dimension2knownRank["5x7x8"] = 206;
    dimension2knownRank["5x7x9"] = 231;
    dimension2knownRank["6x6x10"] = 252;
    dimension2knownRank["7x7x7"] = 250;
    dimension2knownRank["7x7x8"] = 279;
    dimension2knownRank["7x7x9"] = 316;
    dimension2knownRank["7x8x8"] = 310;
    dimension2knownRank["8x8x8"] = 343;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initializeBestBinaryRanks() {
    dimension2knownRank["2x4x5"] = 33;
    dimension2knownRank["2x4x9"] = 59;
    dimension2knownRank["2x4x10"] = 65;
    dimension2knownRank["2x4x13"] = 84;
    dimension2knownRank["2x5x10"] = 80;
    dimension2knownRank["2x7x9"] = 100;
    dimension2knownRank["3x3x6"] = 42;
    dimension2knownRank["3x3x10"] = 71;
    dimension2knownRank["3x3x11"] = 78;
    dimension2knownRank["3x3x12"] = 84;
    dimension2knownRank["3x3x13"] = 91;
    dimension2knownRank["3x3x14"] = 98;
    dimension2knownRank["3x3x15"] = 105;
    dimension2knownRank["3x3x16"] = 112;
    dimension2knownRank["3x4x7"] = 64;
    dimension2knownRank["3x4x13"] = 118;
    dimension2knownRank["3x4x14"] = 127;
    dimension2knownRank["3x4x15"] = 137;
    dimension2knownRank["3x6x6"] = 83;
    dimension2knownRank["3x6x7"] = 96;
    dimension2knownRank["3x6x9"] = 122;
    dimension2knownRank["3x6x10"] = 136;
    dimension2knownRank["3x7x8"] = 128;
    dimension2knownRank["3x7x9"] = 143;
    dimension2knownRank["4x4x4"] = 47;
    dimension2knownRank["4x4x5"] = 60;
    dimension2knownRank["4x4x8"] = 94;
    dimension2knownRank["4x4x9"] = 107;
    dimension2knownRank["4x4x11"] = 132;
    dimension2knownRank["4x4x12"] = 141;
    dimension2knownRank["4x4x13"] = 154;
    dimension2knownRank["4x4x14"] = 167;
    dimension2knownRank["4x4x15"] = 179;
    dimension2knownRank["4x4x16"] = 188;
    dimension2knownRank["4x5x5"] = 73;
    dimension2knownRank["4x5x6"] = 89;
    dimension2knownRank["4x5x9"] = 133;
    dimension2knownRank["4x5x10"] = 146;
    dimension2knownRank["4x5x11"] = 162;
    dimension2knownRank["4x5x12"] = 177;
    dimension2knownRank["4x7x9"] = 187;
    dimension2knownRank["5x5x9"] = 166;
    dimension2knownRank["5x5x10"] = 183;
    dimension2knownRank["5x5x11"] = 200;
    dimension2knownRank["5x5x12"] = 217;
    dimension2knownRank["6x6x10"] = 252;
    dimension2knownRank["7x7x7"] = 248;
    dimension2knownRank["7x7x8"] = 273;
    dimension2knownRank["7x7x9"] = 313;
    dimension2knownRank["7x8x8"] = 302;
    dimension2knownRank["8x8x8"] = 329;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::run() {
    initialize();

    auto startTime = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    for (size_t iteration = 0; 1; iteration++) {
        runIteration();
        updateBest(iteration);
        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);

        report(iteration + 1, startTime, elapsedTimes);

        t1 = std::chrono::high_resolution_clock::now();

        if (resizeProbability > 0) {
            resizeIteration();
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
void MetaFlipGraph<Scheme>::runIteration()  {
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < count; i++)
        randomWalk(schemes[i], schemesBest[i], flips[i], iterations[i], plusIterations[i], bestRanks[i], generators[omp_get_thread_num()]);
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::resizeIteration() {
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < count; i++)
        resize(schemes[i], flips[i], iterations[i], plusIterations[i], generators[omp_get_thread_num()]);    
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

        if (!schemesBest[top].validate()) {
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

        if (!schemes[pair.second].validate()) {
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

    std::cout << "+-----------------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| " << std::left;
    std::cout << "ring: " << std::setw(21) << schemes[0].getRing() << "   ";
    std::cout << "count: " << std::setw(20) << count << "   ";
    std::cout << std::right << std::setw(39) << ("iteration: " + std::to_string(iteration));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "seed: " << std::setw(21) << seed << "   ";
    std::cout << "threads: " << std::setw(18) << threads << "   ";
    std::cout << std::right << std::setw(39) << ("elapsed: " + prettyTime(elapsed));
    std::cout << " |" << std::endl;

    bool improved = false;
    for (const auto &pair : dimension2bestRank) {
        std::string dimension = pair.first;
        int rank = pair.second;

        auto known = dimension2knownRank.find(dimension);
        if (known == dimension2knownRank.end() || rank >= known->second)
            continue;

        if (!improved) {
            std::cout << "+=====================================================================================================+" << std::endl;
            std::cout << "| Improvements:                                                                                       |" << std::endl;
            std::cout << "| +-----------+------------+---------------+                                                          |" << std::endl;
            std::cout << "| | dimension | known rank | improved rank |                                                          |" << std::endl;
            std::cout << "| +-----------+------------+---------------+                                                          |" << std::endl;
            improved = true;
        }

        std::cout << "| | " << std::setw(9) << dimension << " | " << std::setw(10) << known->second << " | " << std::setw(13) << rank << " |                                                          |" << std::endl;
    }

    if (improved) {
        std::cout << "| +-----------+------------+---------------+                                                          |" << std::endl;
        std::cout << "|                                                                                                     |" << std::endl;
    }

    std::cout << "+=====================================================================================================+" << std::endl;
    std::cout << "| runner |   scheme size   | scheme rank |   naive    |            |        flips        |    plus    |" << std::endl;
    std::cout << "|   id   | sorted |  real  | best | curr | complexity | iterations |  count  | available | iterations |" << std::endl;
    std::cout << "+--------+--------+--------+------+------+------------+------------+---------+-----------+------------+" << std::endl;
    std::cout << std::right;

    for (const auto &dimension : dimensions) {
        const std::vector<int> &indices = dimension2indices.at(dimension);

        for (size_t i = 0; i < topCount && i < indices.size(); i++) {
            int runner = indices[i];

            std::cout << "| ";
            std::cout << std::setw(6) << runner << " | ";
            std::cout << std::setw(6) << dimension << " | ";
            std::cout << std::setw(6) << schemes[runner].getDimension() << " | ";
            std::cout << std::setw(4) << bestRanks[runner] << " | ";
            std::cout << std::setw(4) << schemes[runner].getRank() << " | ";
            std::cout << std::setw(10) << schemes[runner].getComplexity() << " | ";
            std::cout << std::setw(10) << prettyInt(iterations[runner]) << " | ";
            std::cout << std::setw(7) << prettyInt(flips[runner]) << " | ";
            std::cout << std::setw(9) << schemes[runner].getAvailableFlips() << " | ";
            std::cout << std::setw(10) << prettyInt(plusIterations[runner]) << " |";
            std::cout << std::endl;
        }

        std::cout << "+--------+--------+--------+------+------+------------+------------+---------+-----------+------------+" << std::endl;
    }

    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator) {
    plusIterations = plusDistribution(generator);

    for (size_t iteration = 0; iteration < flipIterations; iteration++) {
        int prevRank = scheme.getRank();

        if (!scheme.tryFlip(generator)) {
            if (scheme.tryExpand(generator))
                flipsCount = 0;

            continue;
        }

        if (reduceProbability && uniform(generator) < reduceProbability && scheme.tryReduce())
            flipsCount = 0;

        if (sandwichingProbability && uniform(generator) < sandwichingProbability)
            scheme.trySandwiching(generator);

        int rank = scheme.getRank();
        if (rank < prevRank) {
            flipsCount = 0;
            iterationsCount = 0;
        }

        flipsCount++;
        iterationsCount++;

        if (rank < bestRank) {
            schemeBest.copy(scheme);
            bestRank = rank;
        }

        if (flipsCount >= plusIterations && rank < bestRank + plusDiff && scheme.tryExpand(generator))
            flipsCount = 0;

        if (iterationsCount >= resetIterations) {
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
void MetaFlipGraph<Scheme>::resize(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::mt19937 &generator) {
    if (uniform(generator) > resizeProbability)
        return;

    if (uniform(generator) < 0.5)
        scheme.swapSizes(generator);

    int index = generator() % schemesBest.size();
    int minN = 3;
    int maxN = 16;
    int maxRank = 345;
    bool resized = true;

    if (!scheme.tryMerge(schemesBest[index], generator, maxN, maxRank)) {
        double p = uniform(generator);

        if (p < 0.5) {
            resized = scheme.tryProject(generator, minN);
        }
        else {
            resized = scheme.tryExtend(generator, maxN, maxRank);
        }
    }

    if (!resized)
        return;

    flipsCount = 0;
    iterationsCount = 0;
    plusIterations = plusDistribution(generator);
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