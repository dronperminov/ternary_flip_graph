#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>

#include "utils.h"
#include "known_ranks.h"
#include "entities/schemes_rank_pool.hpp"
#include "parameters/flip_parameters.h"
#include "parameters/meta_pool_parameters.h"
#include "parameters/meta_parameters.h"

template <typename Scheme>
class MetaFlipGraphPool {
    int count;
    std::string outputPath;
    int threads;
    FlipParameters flipParameters;
    MetaPoolParameters poolParameters;
    MetaParameters metaParameters;
    int seed;
    std::string format;

    std::vector<std::string> dimensions;
    std::unordered_map<std::string, SchemesRankPool<Scheme>> dimension2pools;
    std::unordered_map<std::string, int> dimension2knownRank;
    std::unordered_map<std::string, double> dimension2priority;

    std::vector<Scheme> schemes;
    std::vector<int> ranks;
    std::vector<size_t> flips;
    std::vector<size_t> iterations;
    std::vector<size_t> plusIterations;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
    std::uniform_int_distribution<size_t> plusDistribution;
public:
    MetaFlipGraphPool(int count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const MetaPoolParameters &poolParameters, const MetaParameters &metaParameters, int seed, const std::string &format);

    bool initializeNaive(int n1, int n2, int n3);
    bool initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness);
    void initializeKnownRanks(const std::string &ring);

    void run();
private:
    bool resume();
    void runIteration();

    void randomWalk(Scheme &scheme, size_t &flipsCount, int &runnerRank, size_t &iterationsCount, size_t &plusIterations, std::vector<Scheme> &pool, std::mt19937 &generator);
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void showImprovements() const;

    void readPriorities();
    void selectRunner(Scheme &scheme, std::mt19937 &generator);
    std::string selectDimension(std::mt19937 &generator);
    void addScheme(const Scheme &scheme, bool save);
    void metaScheme(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator);

    void tryExtend(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator);
    void tryProject(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator);
    void tryProduct(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator);
    void tryMerge(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator);

    bool compareDimension(const std::string &d1, const std::string &d2) const;
    bool canExtend(const std::string &ring, const std::string &dimension, int rank) const;
};

template <typename Scheme>
MetaFlipGraphPool<Scheme>::MetaFlipGraphPool(int count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const MetaPoolParameters &poolParameters, const MetaParameters &metaParameters, int seed, const std::string &format) : uniform(0.0, 1.0), plusDistribution(flipParameters.minPlusIterations, flipParameters.maxPlusIterations) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, count);

    this->flipParameters = flipParameters;
    this->poolParameters = poolParameters;
    this->metaParameters = metaParameters;

    this->seed = seed;
    this->format = format;

    generators = initRandomGenerators(seed, threads);

    schemes.resize(count);
    flips.resize(count);
    ranks.resize(count);
    iterations.resize(count);
    plusIterations.resize(count);
}

template <typename Scheme>
bool MetaFlipGraphPool<Scheme>::initializeNaive(int n1, int n2, int n3) {
    std::cout << "Start initializing with naive " << n1 << "x" << n2 << "x" << n3 << " scheme" << std::endl;

    Scheme scheme;
    if (!scheme.initializeNaive(n1, n2, n3))
        return false;

    addScheme(scheme, false);
    return true;
}

template <typename Scheme>
bool MetaFlipGraphPool<Scheme>::initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    int schemesCount = 1;
    if (multiple)
        f >> schemesCount;

    std::cout << "Start reading " << schemesCount << " schemes from \"" << path << "\"" << std::endl;
    bool valid = true;

    for (int i = 0; i < schemesCount; i++) {
        Scheme scheme;
        if (!scheme.read(f, checkCorrectness)) {
            valid = false;
            break;
        }

        addScheme(scheme, false);
    }

    f.close();
    return valid;
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::initializeKnownRanks(const std::string &ring) {
    dimension2knownRank = KNOWN_RANKS.at(ring);
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::run() {
    if (poolParameters.resume && !resume())
        return;

    iterations.assign(count, 0);

    auto startTime = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    for (size_t iteration = 0; 1; iteration++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        runIteration();
        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);

        report(iteration + 1, startTime, elapsedTimes);
    }
}

template <typename Scheme>
bool MetaFlipGraphPool<Scheme>::resume() {
    if (!std::filesystem::exists(outputPath)) {
        std::cout << "Unable to resume: output directory does not exists (" << outputPath << ")" << std::endl;
        return true;
    }

    std::vector<std::string> paths;
    for (auto it = std::filesystem::recursive_directory_iterator(outputPath); it != std::filesystem::recursive_directory_iterator(); it++)
        if (it->is_regular_file())
            paths.push_back(it->path().string());

    std::cout << "Start adding " << paths.size() << " schemes from " << outputPath << std::endl;
    std::vector<std::vector<Scheme>> pool(threads);

    bool valid = true;

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < paths.size(); i++) {
        Scheme scheme;
        if (!scheme.read(paths[i], false))
            valid = false;

        int thread = omp_get_thread_num();
        pool[thread].emplace_back(scheme);

        if (thread == 0 && pool[thread].size() % 1000 == 0)
            std::cout << "Thread 1 read " << pool[thread].size() << " schemes" << std::endl;
    }

    if (valid) {
        for (int i = 0; i < threads; i++)
            for (const Scheme &scheme : pool[i])
                addScheme(scheme, false);

        std::cout << "All schemes have been read" << std::endl;
    }

    return valid;
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::runIteration() {
    std::vector<std::vector<Scheme>> pool(threads);
    readPriorities();

    for (auto& pair : dimension2pools)
        pair.second.resetDiff();

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++) {
        int thread = omp_get_thread_num();
        randomWalk(schemes[i], flips[i], ranks[i], iterations[i], plusIterations[i], pool[thread], generators[thread]);
    }

    for (int i = 0; i < threads; i++)
        for (const Scheme &scheme : pool[i])
            addScheme(scheme, true);
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::randomWalk(Scheme &scheme, size_t &flipsCount, int &runnerRank, size_t &iterationsCount, size_t &plusIterations, std::vector<Scheme> &pool, std::mt19937 &generator) {
    for (size_t iteration = 0; iteration < flipParameters.flipIterations; iteration++) {
        if (iterationsCount == 0 || iterationsCount >= flipParameters.resetIterations) {
            selectRunner(scheme, generator);
            flipsCount = 0;
            runnerRank = scheme.getRank();
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

        if (scheme.getRank() == runnerRank - 1 && (!poolParameters.liftOnly || scheme.canLift())) {
            Scheme poolScheme;
            poolScheme.copy(scheme);
            pool.emplace_back(poolScheme);
            metaScheme(scheme, pool, generator);
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

        if (flipsCount >= plusIterations && rank < runnerRank + flipParameters.plusDiff && scheme.tryExpand(generator))
            flipsCount = 0;
    }

    if (scheme.getRank() > runnerRank)
        return;

    if (iterationsCount > 0 && uniform(generator) < poolParameters.alternativesProbability && (!poolParameters.liftOnly || scheme.canLift())) {
        Scheme poolScheme;
        poolScheme.copy(scheme);
        pool.emplace_back(poolScheme);
        metaScheme(scheme, pool, generator);
    }
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    double meanFillRatio = 0;
    size_t rediscovered = 0;
    size_t total = 0;
    size_t diff = 0;

    std::cout << "+------------------------------------------------------------------------+" << std::endl << std::left;
    std::cout << "| seed: " << std::setw(64) << seed << " |" << std::endl;
    std::cout << "| runners: " << std::setw(61) << count << " |" << std::endl;
    std::cout << "| threads: " << std::setw(61) << threads << " |" << std::endl;
    std::cout << "| iteration: " << std::setw(59) << iteration << " |" << std::endl;
    std::cout << "| elapsed time: " << std::setw(56) << prettyTime(elapsed) << " |" << std::endl;
    std::cout << "+-----------+------+-----------------+---------------+-------------------+" << std::endl;
    std::cout << "|           |      |                 |  complexity   |  available flips  |" << std::endl;
    std::cout << "| dimension | rank |   total count   |  min  |  max  |   min   |   max   |" << std::endl;

    for (const std::string &dimension : dimensions) {
        const SchemesRankPool<Scheme> &pool = dimension2pools.at(dimension);
        int knownRank = dimension2knownRank.at(dimension);
        diff += pool.print(knownRank);
        total += pool.size();

        meanFillRatio += pool.minFillRatio();
        if (pool.minRank() == knownRank)
            rediscovered++;
    }

    std::cout << "+-----------+------+-----------------+---------------+-------------------+" << std::endl;
    showImprovements();
    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << "- mean fill ratio: " << std::setprecision(3) << (meanFillRatio / dimensions.size()) << std::endl;
    std::cout << "- total schemes: " << total;
    if (diff)
        std::cout << " (+" << diff << ")";
    std::cout << std::endl;
    std::cout << "- rediscovered: " << rediscovered << " / " << dimensions.size() << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::showImprovements() const {
    bool showed = false;

    for (const std::string &dimension : dimensions) {
        int known = dimension2knownRank.at(dimension);
        int curr = dimension2pools.at(dimension).minRank();
        if (curr >= known)
            continue;

        if (!showed) {
            showed = true;
            std::cout << "| improvements:                                                          |" << std::endl;
        }

        std::stringstream ss;
        ss << dimension << ": " << curr << " (" << known << ")";
        std::cout << "| " << std::setw(70) << ss.str() << " |" << std::endl;
    }

    if (showed)
        std::cout << "+------------------------------------------------------------------------+" << std::endl << std::left;
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::readPriorities() {
    dimension2priority.clear();

    if (!std::filesystem::exists(poolParameters.prioritiesPath))
        return;

    std::ifstream f(poolParameters.prioritiesPath);
    if (!f)
        return;

    std::string line;

    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::string dimension;
        double priority;
        std::stringstream ss(line);

        if (!(ss >> dimension >> priority)) {
            std::cout << "Invalid priority line: " << line << std::endl;
            continue;
        }

        if (dimension2knownRank.find(dimension) == dimension2knownRank.end()) {
            std::cout << "Invalid priority dimension: " << dimension << std::endl;
            continue;
        }

        dimension2priority[dimension] = priority;
        std::cout << "Set priority " << priority << " for " << dimension << std::endl;
    }

    f.close();
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::selectRunner(Scheme &scheme, std::mt19937 &generator) {
    std::string dimension = selectDimension(generator);
    SchemesRankPool<Scheme> &pools = dimension2pools.at(dimension);
    pools.copyRandom(scheme, generator, poolParameters.selectRankScale);
}

template <typename Scheme>
std::string MetaFlipGraphPool<Scheme>::selectDimension(std::mt19937 &generator) {
    std::vector<double> weights(dimensions.size());
    double total = 0;

    for (size_t i = 0; i < dimensions.size(); i++) {
        const SchemesRankPool<Scheme> &pool = dimension2pools.at(dimensions[i]);
        int rank = std::min(dimension2knownRank.at(dimensions[i]), pool.minRank());

        if (dimension2priority.find(dimensions[i]) != dimension2priority.end())
            weights[i] = dimension2priority.at(dimensions[i]);
        else
            weights[i] = 1.0 - pool.fillRatio(rank) + 0.1 / dimensions.size();

        total += weights[i];
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double randomWeight = uniform(generator) * total;
    double sum = 0;

    for (size_t i = 0; i < dimensions.size(); i++) {
        sum += weights[i];

        if (randomWeight <= sum)
            return dimensions[i];
    }

    return dimensions.back();
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::addScheme(const Scheme &scheme, bool save) {
    std::string dimension = scheme.getDimension();

    if (dimension2pools.find(dimension) == dimension2pools.end()) {
        dimensions.push_back(dimension);
        std::sort(dimensions.begin(), dimensions.end(), [&](const std::string &d1, const std::string &d2) { return compareDimension(d1, d2); });
        dimension2pools.emplace(dimension, SchemesRankPool<Scheme>(dimension, poolParameters.size, poolParameters.uniqueOnly, outputPath + "/" + dimension, format));
    }

    dimension2pools.at(dimension).add(scheme, save);
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::metaScheme(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator) {
    std::string dimension = scheme.getDimension();
    int rank = scheme.getRank();

    const auto& pool = dimension2pools.at(dimension);

    int minRank = std::min(pool.minRank(), dimension2knownRank.at(dimension));

    if (rank > minRank + metaParameters.maxRankDiff || !canExtend(scheme.getRing(), dimension, rank) || uniform(generator) >= metaParameters.probability * pow(poolParameters.metaRankScale, rank - minRank))
        return;

    if (uniform(generator) < poolParameters.extendProbability)
        tryExtend(scheme, schemesPool, generator);

    if (uniform(generator) < poolParameters.projectProbability)
        tryProject(scheme, schemesPool, generator);

    if (uniform(generator) < poolParameters.productProbability)
        tryProduct(scheme, schemesPool, generator);

    if (uniform(generator) < poolParameters.mergeProbability)
        tryMerge(scheme, schemesPool, generator);
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::tryExtend(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator) {
    for (int i = 0; i < 3; i++) {
        if (!scheme.isValidExtension(i, metaParameters.maxDimension, metaParameters.maxRank))
            continue;

        Scheme poolScheme;
        poolScheme.copy(scheme);
        poolScheme.extend(i);
        poolScheme.fixSizes();

        if (poolScheme.getRank() <= dimension2knownRank.at(poolScheme.getDimension()) + poolParameters.extendMaxDiff && (!poolParameters.liftOnly || poolScheme.canLift()))
            schemesPool.emplace_back(poolScheme);
    }
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::tryProject(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator) {
    int n1 = scheme.getDimension(0);
    int n2 = scheme.getDimension(1);
    int n3 = scheme.getDimension(2);

    if (n1 < poolParameters.projectMinN1 || n2 < poolParameters.projectMinN2 || n3 < poolParameters.projectMinN3)
        return;

    for (int i = 0; i < 3; i++) {
        if (!scheme.isValidProject(i, metaParameters.minDimension))
            continue;

        for (int j = 0; j < scheme.getDimension(i); j++) {
            Scheme poolScheme;
            poolScheme.copy(scheme);
            poolScheme.project(i, j);
            poolScheme.fixSizes();

            if (poolScheme.getRank() <= dimension2knownRank.at(poolScheme.getDimension()) + poolParameters.projectMaxDiff && (!poolParameters.liftOnly || poolScheme.canLift()))
                schemesPool.emplace_back(poolScheme);
        }
    }
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::tryProduct(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator) {
    for (int n1 = 2; scheme.getDimension(0) * n1 <= metaParameters.maxDimension; n1++) {
        for (int n2 = 2; scheme.getDimension(1) * n2 <= metaParameters.maxDimension; n2++) {
            for (int n3 = 2; scheme.getDimension(2) * n3 <= metaParameters.maxDimension; n3++) {
                std::string dimension = getDimension(n1, n2, n3);
                if (dimension2pools.find(dimension) == dimension2pools.end())
                    continue;

                Scheme scheme2;
                dimension2pools.at(dimension).copyRandomMinRank(scheme2, generator);
                scheme2.setSizes(n1, n2, n3);

                if (!scheme.isValidProduct(scheme2, metaParameters.maxDimension, metaParameters.maxRank))
                    continue;

                Scheme poolScheme;
                poolScheme.copy(scheme);
                poolScheme.product(scheme2);

                if (poolScheme.getRank() <= dimension2knownRank.at(poolScheme.getDimension()) + poolParameters.productMaxDiff && (!poolParameters.liftOnly || poolScheme.canLift()))
                    schemesPool.emplace_back(poolScheme);
            }
        }
    }
}

template <typename Scheme>
void MetaFlipGraphPool<Scheme>::tryMerge(const Scheme &scheme, std::vector<Scheme> &schemesPool, std::mt19937 &generator) {
    for (int i = 0; i < 3; i++) {
        for (int ni = 2; ni <= scheme.getDimension(i) && ni + scheme.getDimension(i) <= metaParameters.maxDimension; ni++) {
            int n[3] = {scheme.getDimension(0), scheme.getDimension(1), scheme.getDimension(2)};
            n[i] = ni;

            std::string dimension = getDimension(n[0], n[1], n[2], true);
            if (dimension2pools.find(dimension) == dimension2pools.end())
                continue;

            const auto& pool = dimension2pools.at(dimension);

            Scheme scheme2;
            pool.copyRandomMinRank(scheme2, generator);
            scheme2.setSizes(n[0], n[1], n[2]);

            if (!scheme.isValidMerge(i, scheme2, metaParameters.maxDimension, metaParameters.maxRank))
                continue;

            Scheme poolScheme;
            poolScheme.copy(scheme);
            poolScheme.merge(scheme2, i);
            poolScheme.fixSizes();

            if (poolScheme.getRank() <= dimension2knownRank.at(poolScheme.getDimension()) + poolParameters.mergeMaxDiff && (!poolParameters.liftOnly || poolScheme.canLift()))
                schemesPool.emplace_back(poolScheme);
        }
    }
}

template <typename Scheme>
bool MetaFlipGraphPool<Scheme>::compareDimension(const std::string &d1, const std::string &d2) const {
    std::stringstream ss1(d1);
    std::stringstream ss2(d2);
    int n1, n2;
    char c;

    for (int i = 0; i < 3; i++) {
        ss1 >> n1 >> c;
        ss2 >> n2 >> c;
        if (n1 != n2)
            return n1 < n2;
    }

    return false;
}

template <typename Scheme>
bool MetaFlipGraphPool<Scheme>::canExtend(const std::string &ring, const std::string &dimension, int rank) const {
    if (ring != "Z2")
        return true;

    std::unordered_set<std::string> ignored = {
        "3x3x8", "3x3x12", "3x3x13", "3x3x14", "3x3x15", "3x3x16", "4x4x4", "4x4x5", "4x5x5"
    };

    return ignored.find(dimension) == ignored.end() || rank >= dimension2knownRank.at(dimension);
}
