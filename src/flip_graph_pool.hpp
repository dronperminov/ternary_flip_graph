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
#include "entities/flip_parameters.h"
#include "entities/pool_parameters.h"

template <typename Scheme>
class FlipGraphPool {
    int count;
    std::string outputPath;
    int threads;
    FlipParameters flipParameters;
    PoolParameters poolParameters;

    int seed;
    int topCount;
    std::string format;

    std::vector<Scheme> initialPool;
    std::vector<Scheme> pool;
    std::vector<Scheme> schemes;
    std::vector<size_t> flips;
    std::vector<size_t> iterations;
    std::vector<size_t> plusIterations;
    std::vector<int> indices;
    int poolRank;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
    std::uniform_int_distribution<size_t> plusDistribution;
public:
    FlipGraphPool(int count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const PoolParameters &poolParameters, int seed, int topCount, const std::string &format);

    bool initializeNaive(int n1, int n2, int n3);
    bool initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness);

    void run(int targetRank);
private:
    void initIteration();
    void runIteration();
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes);
    void updatePool();

    void randomWalk(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::vector<Scheme> &pool, std::mt19937 &generator);

    std::string getSavePath(const Scheme &scheme, int version, const std::string path) const;
    void saveScheme(const Scheme &scheme, const std::string &path) const;
};

template <typename Scheme>
FlipGraphPool<Scheme>::FlipGraphPool(int count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const PoolParameters &poolParameters, int seed, int topCount, const std::string &format) : uniform(0.0, 1.0), plusDistribution(flipParameters.minPlusIterations, flipParameters.maxPlusIterations) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, count);
    this->flipParameters = flipParameters;
    this->poolParameters = poolParameters;
    this->seed = seed;
    this->topCount = std::min(topCount, count);
    this->format = format;

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

    schemes.resize(count);
    flips.resize(count);
    iterations.resize(count);
    plusIterations.resize(count);
    indices.resize(count);

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        indices[i] = i;
}

template <typename Scheme>
bool FlipGraphPool<Scheme>::initializeNaive(int n1, int n2, int n3) {
    std::cout << "Start initializing pool with naive " << n1 << "x" << n2 << "x" << n3 << " scheme" << std::endl;

    Scheme scheme;
    if (!scheme.initializeNaive(n1, n2, n3))
        return false;

    initialPool.emplace_back(scheme);
    poolRank = n1 * n2 * n3;

    std::cout << "Initial pool rank set to " << poolRank << std::endl;
    return true;
}

template <typename Scheme>
bool FlipGraphPool<Scheme>::initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    int schemesCount = 1;
    if (multiple)
        f >> schemesCount;

    int readCount = std::min(schemesCount, (int)poolParameters.size);
    std::cout << "Start reading " << readCount << " / " << schemesCount << " schemes from \"" << path << "\" as initial pool" << std::endl;

    initialPool.resize(readCount);

    bool valid = true;

    for (int i = 0; i < readCount; i++) {
        if (!initialPool[i].read(f, checkCorrectness)) {
            valid = false;
            break;
        }

        int schemeRank = initialPool[i].getRank();

        if (i == 0) {
            poolRank = schemeRank;
            std::cout << "Initial pool rank set to " << poolRank << std::endl;
            continue;
        }

        if (poolRank != schemeRank) {
            std::cout << "Scheme " << (i + 1) << " rank is different (" << schemeRank << " != " << poolRank << ")" << std::endl;
            valid = false;
            break;
        }
    }

    f.close();
    return valid;
}

template <typename Scheme>
void FlipGraphPool<Scheme>::run(int targetRank) {
    auto startTime = std::chrono::high_resolution_clock::now();

    while (poolRank > targetRank) {
        initIteration();
        std::vector<double> elapsedTimes;

        for (size_t iteration = 0; pool.size() < poolParameters.size; iteration++) {
            auto t1 = std::chrono::high_resolution_clock::now();
            runIteration();
            auto t2 = std::chrono::high_resolution_clock::now();
            elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);

            report(iteration + 1, startTime, elapsedTimes);
        }

        updatePool();
    }
}

template <typename Scheme>
void FlipGraphPool<Scheme>::initIteration() {
    std::cout << "Run next iteration with pool rank " << poolRank << " (" << initialPool.size() << " schemes)" << std::endl;
    poolRank--;

    pool.clear();
    iterations.assign(count, 0);

    std::string poolPath = outputPath + "/rank" + std::to_string(poolRank);
    if (!makeDirectory(poolPath))
        exit(-1);
}

template <typename Scheme>
void FlipGraphPool<Scheme>::runIteration() {
    std::vector<std::vector<Scheme>> poolIteration(threads);

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++) {
        int thread = omp_get_thread_num();
        randomWalk(schemes[i], flips[i], iterations[i], plusIterations[i], poolIteration[thread], generators[thread]);
    }

    std::string path = outputPath + "/rank" + std::to_string(poolRank);
    for (int i = 0; i < threads; i++) {
        for (const Scheme &scheme : poolIteration[i]) {
            pool.emplace_back(scheme);
            std::string schemePath = getSavePath(scheme, pool.size(), path);
            saveScheme(scheme, schemePath);
        }
    }
}

template <typename Scheme>
void FlipGraphPool<Scheme>::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    std::cout << std::left;
    std::cout << "+------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| ";
    std::cout << "dimension: " << std::setw(14) << initialPool[0].getDimension() << "   ";
    std::cout << "seed: " << std::setw(20) << seed << "   ";
    std::cout << std::right << std::setw(19) << ("iteration: " + std::to_string(iteration));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "threads: " << std::setw(16) << threads << "   ";
    std::cout << "flip iters: " << std::setw(14) << prettyInt(flipParameters.flipIterations) << "   ";
    std::cout << std::right << std::setw(19) << ("elapsed: " + prettyTime(elapsed));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "runners: " << std::setw(16) << count << "   ";
    std::cout << "reset iters: " << std::setw(13) << prettyInt(flipParameters.resetIterations) << "   ";
    std::cout << std::right << std::setw(19) << ("pool size: " + std::to_string(pool.size()));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "ring: " << std::setw(19) << initialPool[0].getRing() << "   ";
    std::cout << "plus diff: " << std::setw(15) << flipParameters.plusDiff << "   ";
    std::cout << std::right << std::setw(19) << ("pool rank: " + std::to_string(poolRank));
    std::cout << " |" << std::endl;

    std::partial_sort(indices.begin(), indices.begin() + topCount, indices.end(), [this](int index1, int index2) {
        int rank1 = schemes[index1].getRank();
        int rank2 = schemes[index2].getRank();
        if (rank1 != rank2)
            return rank1 < rank2;

        int complexity1 = schemes[index1].getComplexity();
        int complexity2 = schemes[index2].getComplexity();
        return complexity1 < complexity2;
    });

    std::cout << "+==============================================================================+" << std::endl;
    std::cout << "| runner | runner |   naive    |            |        flips        |    plus    |" << std::endl;
    std::cout << "|   id   |  rank  | complexity | iterations |  count  | available | iterations |" << std::endl;
    std::cout << "+--------+--------+------------+------------+---------+-----------+------------+" << std::endl;
    std::cout << std::right;

    for (int i = 0; i < topCount; i++) {
        int runner = indices[i];
        std::cout << "| ";
        std::cout << std::setw(6) << (runner + 1) << " | ";
        std::cout << std::setw(6) << schemes[runner].getRank() << " | ";
        std::cout << std::setw(10) << schemes[runner].getComplexity() << " | ";
        std::cout << std::setw(10) << prettyInt(iterations[runner]) << " | ";
        std::cout << std::setw(7) << prettyInt(flips[runner]) << " | ";
        std::cout << std::setw(9) << schemes[runner].getAvailableFlips() << " | ";
        std::cout << std::setw(10) << prettyInt(plusIterations[runner]) << " |";
        std::cout << std::endl;
    }

    std::cout << "+--------+--------+------------+------------+---------+-----------+------------+" << std::endl;
    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
void FlipGraphPool<Scheme>::updatePool() {
    initialPool.resize(pool.size());

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < pool.size(); i++)
        initialPool[i].copy(pool[i]);
}

template <typename Scheme>
void FlipGraphPool<Scheme>::randomWalk(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::vector<Scheme> &pool, std::mt19937 &generator) {
    size_t size = this->pool.size();

    for (size_t iteration = 0; iteration < flipParameters.flipIterations; iteration++) {
        if (iterationsCount == 0 || iterationsCount >= flipParameters.resetIterations) {
            scheme.copy(initialPool[generator() % initialPool.size()]);
            flipsCount = 0;
            iterationsCount = 0;
            plusIterations = plusDistribution(generator);
        }

        iterationsCount++;
        int prevRank = scheme.getRank();

        if (!scheme.tryFlip(generator)) {
            if (scheme.tryExpand(generator))
                flipsCount = 0;

            continue;
        }

        if (scheme.getRank() == poolRank) {
            Scheme poolScheme;
            poolScheme.copy(scheme);
            pool.emplace_back(poolScheme);

            if (size + pool.size() >= poolParameters.size)
                break;

            iterationsCount = 0;
            continue;
        }

        flipsCount++;

        if (flipParameters.reduceProbability && uniform(generator) < flipParameters.reduceProbability && scheme.tryReduce())
            flipsCount = 0;

        if (flipParameters.sandwichingProbability && uniform(generator) < flipParameters.sandwichingProbability)
            scheme.trySandwiching(generator);

        int rank = scheme.getRank();
        if (rank < prevRank)
            flipsCount = 0;

        if (flipsCount >= plusIterations && rank < poolRank + 1 + flipParameters.plusDiff && scheme.tryExpand(generator))
            flipsCount = 0;
    }
}

template <typename Scheme>
std::string FlipGraphPool<Scheme>::getSavePath(const Scheme &scheme, int version, const std::string path) const {
    std::stringstream ss;
    ss << path << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_version" << std::setfill('0') << std::setw(5) << version;
    ss << "_c" << scheme.getComplexity();
    ss << "_" << scheme.getRing();
    return ss.str();
}

template <typename Scheme>
void FlipGraphPool<Scheme>::saveScheme(const Scheme &scheme, const std::string &path) const {
    if (format == "json") {
        scheme.saveJson(path + ".json");
    }
    else if (format == "txt") {
        scheme.saveTxt(path + ".txt");
    }
}
