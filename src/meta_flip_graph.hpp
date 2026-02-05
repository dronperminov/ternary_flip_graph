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
#include "entities/flip_parameters.h"
#include "entities/meta_parameters.h"

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
    std::vector<std::string> dimensions;
public:
    MetaFlipGraph(size_t count, const std::string outputPath, int threads, const FlipParameters &flipParameters, const MetaParameters &metaParameters, int seed, size_t topCount, const std::string &format);

    bool initializeNaive(int n1, int n2, int n3);
    bool initializeFromFile(const std::string &path, bool multiple);

    void initializeKnownRanks(const std::string &ring);
    void run();
private:
    void initialize();
    void initializeKnownRationalRanks();
    void initializeKnownTernaryRanks();
    void initializeKnownBinaryRanks();

    void flipIteration();
    void metaIteration();

    void updateBest(size_t iteration);
    void updateRanks(int iteration, bool save);
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void randomWalk(Scheme &scheme, Scheme &schemeBest, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, int &bestRank, std::mt19937 &generator);
    void meta(Scheme &scheme, size_t &flipsCount, size_t &iterationsCount, size_t &plusIterations, std::mt19937 &generator);

    void updateIndices();
    bool compare(int index1, int index2) const;
    std::string getSavePath(const Scheme &scheme, int iteration, const std::string path) const;
    std::string sortedDimension(const Scheme &scheme) const;

    void saveScheme(const Scheme &scheme, const std::string &path) const;
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

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

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
    bool unknown = dimension2knownRank.empty();

    for (size_t i = 0; i < count && i < schemesCount; i++) {
        std::string dimension = sortedDimension(schemes[i]);
        dimension2improvements[dimension].push_back(Scheme(schemes[i]));

        if (dimension2knownRank.find(dimension) == dimension2knownRank.end() || (unknown && schemes[i].getRank() < dimension2knownRank.at(dimension)))
            dimension2knownRank[dimension] = schemes[i].getRank();
    }

    #pragma omp parallel for num_threads(threads)
    for (size_t i = schemesCount; i < count; i++)
        schemes[i].copy(schemes[i % schemesCount]);

    return true;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initializeKnownRanks(const std::string &ring) {
    if (ring == "Q") {
        initializeKnownRationalRanks();
    }
    else if (ring == "ZT") {
        initializeKnownTernaryRanks();
    }
    else if (ring == "Z2") {
        initializeKnownBinaryRanks();
    }
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
void MetaFlipGraph<Scheme>::initializeKnownRationalRanks() {
    dimension2knownRank = {
        {"2x2x2", 7}, {"2x2x3", 11}, {"2x2x4", 14}, {"2x2x5", 18}, {"2x2x6", 21}, {"2x2x7", 25}, {"2x2x8", 28}, {"2x2x9", 32}, {"2x2x10", 35}, {"2x2x11", 39}, {"2x2x12", 42}, {"2x2x13", 46}, {"2x2x14", 49}, {"2x2x15", 53}, {"2x2x16", 56},
        {"2x3x3", 15}, {"2x3x4", 20}, {"2x3x5", 25}, {"2x3x6", 30}, {"2x3x7", 35}, {"2x3x8", 40}, {"2x3x9", 45}, {"2x3x10", 50}, {"2x3x11", 55}, {"2x3x12", 60}, {"2x3x13", 65}, {"2x3x14", 70}, {"2x3x15", 75}, {"2x3x16", 80},
        {"2x4x4", 26}, {"2x4x5", 32}, {"2x4x6", 39}, {"2x4x7", 45}, {"2x4x8", 51}, {"2x4x9", 58}, {"2x4x10", 64}, {"2x4x11", 71}, {"2x4x12", 77}, {"2x4x13", 83}, {"2x4x14", 90}, {"2x4x15", 96}, {"2x4x16", 102},
        {"2x5x5", 40}, {"2x5x6", 47}, {"2x5x7", 55}, {"2x5x8", 63}, {"2x5x9", 72}, {"2x5x10", 79}, {"2x5x11", 87}, {"2x5x12", 94}, {"2x5x13", 102}, {"2x5x14", 110}, {"2x5x15", 118}, {"2x5x16", 126},
        {"2x6x6", 56}, {"2x6x7", 66}, {"2x6x8", 75}, {"2x6x9", 86}, {"2x6x10", 94}, {"2x6x11", 103}, {"2x6x12", 112}, {"2x6x13", 122}, {"2x6x14", 131}, {"2x6x15", 141}, {"2x6x16", 150},
        {"2x7x7", 76}, {"2x7x8", 88}, {"2x7x9", 99}, {"2x7x10", 110}, {"2x7x11", 121}, {"2x7x12", 131}, {"2x7x13", 142}, {"2x7x14", 152}, {"2x7x15", 164}, {"2x7x16", 175},
        {"2x8x8", 100}, {"2x8x9", 113}, {"2x8x10", 125}, {"2x8x11", 138}, {"2x8x12", 150}, {"2x8x13", 164}, {"2x8x14", 175}, {"2x8x15", 188}, {"2x8x16", 200},
        {"2x9x9", 126}, {"2x9x10", 140}, {"2x9x11", 154}, {"2x9x12", 168}, {"2x9x13", 182}, {"2x9x14", 196},
        {"2x10x10", 155}, {"2x10x11", 171}, {"2x10x12", 186},
        {"2x11x11", 187},
        {"3x3x3", 23}, {"3x3x4", 29}, {"3x3x5", 36}, {"3x3x6", 40}, {"3x3x7", 49}, {"3x3x8", 55}, {"3x3x9", 63}, {"3x3x10", 69}, {"3x3x11", 76}, {"3x3x12", 80}, {"3x3x13", 89}, {"3x3x14", 95}, {"3x3x15", 103}, {"3x3x16", 109},
        {"3x4x4", 38}, {"3x4x5", 47}, {"3x4x6", 54}, {"3x4x7", 63}, {"3x4x8", 73}, {"3x4x9", 83}, {"3x4x10", 92}, {"3x4x11", 101}, {"3x4x12", 108}, {"3x4x13", 117}, {"3x4x14", 126}, {"3x4x15", 136}, {"3x4x16", 146},
        {"3x5x5", 58}, {"3x5x6", 68}, {"3x5x7", 79}, {"3x5x8", 90}, {"3x5x9", 104}, {"3x5x10", 115}, {"3x5x11", 126}, {"3x5x12", 136}, {"3x5x13", 147}, {"3x5x14", 158}, {"3x5x15", 169}, {"3x5x16", 180},
        {"3x6x6", 80}, {"3x6x7", 94}, {"3x6x8", 108}, {"3x6x9", 120}, {"3x6x10", 134}, {"3x6x11", 148}, {"3x6x12", 160}, {"3x6x13", 174}, {"3x6x14", 188}, {"3x6x15", 200}, {"3x6x16", 214},
        {"3x7x7", 111}, {"3x7x8", 126}, {"3x7x9", 142}, {"3x7x10", 157}, {"3x7x11", 173}, {"3x7x12", 188}, {"3x7x13", 205}, {"3x7x14", 220}, {"3x7x15", 236}, {"3x7x16", 251},
        {"3x8x8", 145}, {"3x8x9", 163}, {"3x8x10", 180}, {"3x8x11", 198}, {"3x8x12", 216}, {"3x8x13", 234}, {"3x8x14", 252}, {"3x8x15", 270}, {"3x8x16", 288},
        {"3x9x9", 183}, {"3x9x10", 203}, {"3x9x11", 224}, {"3x9x12", 240}, {"3x9x13", 262}, {"3x9x14", 283},
        {"3x10x10", 226}, {"3x10x11", 249}, {"3x10x12", 268},
        {"3x11x11", 274},
        {"4x4x4", 48}, {"4x4x5", 61}, {"4x4x6", 73}, {"4x4x7", 85}, {"4x4x8", 96}, {"4x4x9", 104}, {"4x4x10", 115}, {"4x4x11", 130}, {"4x4x12", 141}, {"4x4x13", 152}, {"4x4x14", 163}, {"4x4x15", 176}, {"4x4x16", 188},
        {"4x5x5", 76}, {"4x5x6", 90}, {"4x5x7", 104}, {"4x5x8", 118}, {"4x5x9", 132}, {"4x5x10", 146}, {"4x5x11", 160}, {"4x5x12", 175}, {"4x5x13", 192}, {"4x5x14", 207}, {"4x5x15", 221}, {"4x5x16", 236},
        {"4x6x6", 105}, {"4x6x7", 123}, {"4x6x8", 140}, {"4x6x9", 159}, {"4x6x10", 175}, {"4x6x11", 194}, {"4x6x12", 210}, {"4x6x13", 228}, {"4x6x14", 245}, {"4x6x15", 263}, {"4x6x16", 280},
        {"4x7x7", 144}, {"4x7x8", 164}, {"4x7x9", 186}, {"4x7x10", 203}, {"4x7x11", 227}, {"4x7x12", 246}, {"4x7x13", 266}, {"4x7x14", 285}, {"4x7x15", 307}, {"4x7x16", 324},
        {"4x8x8", 182}, {"4x8x9", 206}, {"4x8x10", 224}, {"4x8x11", 252}, {"4x8x12", 272}, {"4x8x13", 297}, {"4x8x14", 315}, {"4x8x15", 339}, {"4x8x16", 357},
        {"4x9x9", 225}, {"4x9x10", 255}, {"4x9x11", 279}, {"4x9x12", 300}, {"4x9x13", 329}, {"4x9x14", 355},
        {"4x10x10", 280}, {"4x10x11", 308}, {"4x10x12", 329},
        {"4x11x11", 340},
        {"5x5x5", 93}, {"5x5x6", 110}, {"5x5x7", 127}, {"5x5x8", 144}, {"5x5x9", 163}, {"5x5x10", 184}, {"5x5x11", 202}, {"5x5x12", 220}, {"5x5x13", 237}, {"5x5x14", 254}, {"5x5x15", 271}, {"5x5x16", 288},
        {"5x6x6", 130}, {"5x6x7", 150}, {"5x6x8", 170}, {"5x6x9", 197}, {"5x6x10", 217}, {"5x6x11", 236}, {"5x6x12", 250}, {"5x6x13", 278}, {"5x6x14", 297}, {"5x6x15", 318}, {"5x6x16", 340},
        {"5x7x7", 176}, {"5x7x8", 204}, {"5x7x9", 229}, {"5x7x10", 254}, {"5x7x11", 277}, {"5x7x12", 296}, {"5x7x13", 325}, {"5x7x14", 349}, {"5x7x15", 375}, {"5x7x16", 398},
        {"5x8x8", 230}, {"5x8x9", 260}, {"5x8x10", 284}, {"5x8x11", 312}, {"5x8x12", 333}, {"5x8x13", 363}, {"5x8x14", 387}, {"5x8x15", 419}, {"5x8x16", 445},
        {"5x9x9", 294}, {"5x9x10", 322}, {"5x9x11", 353}, {"5x9x12", 377}, {"5x9x13", 411}, {"5x9x14", 439},
        {"5x10x10", 352}, {"5x10x11", 386}, {"5x10x12", 413},
        {"5x11x11", 424},
        {"6x6x6", 153}, {"6x6x7", 183}, {"6x6x8", 203}, {"6x6x9", 225}, {"6x6x10", 247}, {"6x6x11", 268}, {"6x6x12", 280}, {"6x6x13", 316}, {"6x6x14", 336}, {"6x6x15", 360}, {"6x6x16", 385},
        {"6x7x7", 212}, {"6x7x8", 238}, {"6x7x9", 268}, {"6x7x10", 296}, {"6x7x11", 318}, {"6x7x12", 336}, {"6x7x13", 372}, {"6x7x14", 399}, {"6x7x15", 430}, {"6x7x16", 457},
        {"6x8x8", 266}, {"6x8x9", 296}, {"6x8x10", 329}, {"6x8x11", 357}, {"6x8x12", 378}, {"6x8x13", 414}, {"6x8x14", 441}, {"6x8x15", 480}, {"6x8x16", 511},
        {"6x9x9", 342}, {"6x9x10", 373}, {"6x9x11", 407}, {"6x9x12", 434}, {"6x9x13", 474}, {"6x9x14", 500},
        {"6x10x10", 406}, {"6x10x11", 446}, {"6x10x12", 476},
        {"6x11x11", 490},
        {"7x7x7", 249}, {"7x7x8", 277}, {"7x7x9", 315}, {"7x7x10", 346}, {"7x7x11", 376}, {"7x7x12", 402}, {"7x7x13", 441}, {"7x7x14", 471}, {"7x7x15", 508}, {"7x7x16", 539},
        {"7x8x8", 306}, {"7x8x9", 350}, {"7x8x10", 385}, {"7x8x11", 423}, {"7x8x12", 454}, {"7x8x13", 496}, {"7x8x14", 529}, {"7x8x15", 571}, {"7x8x16", 603},
        {"7x9x9", 398}, {"7x9x10", 437}, {"7x9x11", 480}, {"7x9x12", 510}, {"7x9x13", 562}, {"7x9x14", 597},
        {"7x10x10", 478}, {"7x10x11", 526}, {"7x10x12", 564},
        {"7x11x11", 577},
        {"8x8x8", 336}, {"8x8x9", 388}, {"8x8x10", 427}, {"8x8x11", 475}, {"8x8x12", 504}, {"8x8x13", 559}, {"8x8x14", 595}, {"8x8x15", 635}, {"8x8x16", 672},
        {"8x9x9", 430}, {"8x9x10", 487}, {"8x9x11", 533}, {"8x9x12", 560}, {"8x9x13", 624}, {"8x9x14", 669},
        {"8x10x10", 532}, {"8x10x11", 588}, {"8x10x12", 630},
        {"8x11x11", 641},
        {"9x9x9", 498}, {"9x9x10", 534}, {"9x9x11", 576}, {"9x9x12", 600}, {"9x9x13", 681}, {"9x9x14", 726},
        {"9x10x10", 600}, {"9x10x11", 651}, {"9x10x12", 684},
        {"9x11x11", 725},
        {"10x10x10", 651}, {"10x10x11", 719}, {"10x10x12", 770},
        {"10x11x11", 793},
        {"11x11x11", 873}
    };

    std::cout << "Initialized known Q ranks" << std::endl;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initializeKnownTernaryRanks() {
    dimension2knownRank = {
        {"2x2x2", 7}, {"2x2x3", 11}, {"2x2x4", 14}, {"2x2x5", 18}, {"2x2x6", 21}, {"2x2x7", 25}, {"2x2x8", 28}, {"2x2x9", 32}, {"2x2x10", 35}, {"2x2x11", 39}, {"2x2x12", 42}, {"2x2x13", 46}, {"2x2x14", 49}, {"2x2x15", 53}, {"2x2x16", 56},
        {"2x3x3", 15}, {"2x3x4", 20}, {"2x3x5", 25}, {"2x3x6", 30}, {"2x3x7", 35}, {"2x3x8", 40}, {"2x3x9", 45}, {"2x3x10", 50}, {"2x3x11", 55}, {"2x3x12", 60}, {"2x3x13", 65}, {"2x3x14", 70}, {"2x3x15", 75}, {"2x3x16", 80},
        {"2x4x4", 26}, {"2x4x5", 33}, {"2x4x6", 39}, {"2x4x7", 45}, {"2x4x8", 51}, {"2x4x9", 59}, {"2x4x10", 65}, {"2x4x11", 71}, {"2x4x12", 77}, {"2x4x13", 84}, {"2x4x14", 90}, {"2x4x15", 96}, {"2x4x16", 102},
        {"2x5x5", 40}, {"2x5x6", 47}, {"2x5x7", 57}, {"2x5x8", 65}, {"2x5x9", 72}, {"2x5x10", 80}, {"2x5x11", 87}, {"2x5x12", 94}, {"2x5x13", 104}, {"2x5x14", 112}, {"2x5x15", 119}, {"2x5x16", 127},
        {"2x6x6", 57}, {"2x6x7", 67}, {"2x6x8", 77}, {"2x6x9", 86}, {"2x6x10", 94}, {"2x6x11", 104}, {"2x6x12", 114}, {"2x6x13", 124}, {"2x6x14", 133}, {"2x6x15", 141}, {"2x6x16", 151},
        {"2x7x7", 77}, {"2x7x8", 88}, {"2x7x9", 102}, {"2x7x10", 112}, {"2x7x11", 122}, {"2x7x12", 133}, {"2x7x13", 144}, {"2x7x14", 154}, {"2x7x15", 165}, {"2x7x16", 176},
        {"2x8x8", 100}, {"2x8x9", 116}, {"2x8x10", 128}, {"2x8x11", 139}, {"2x8x12", 151}, {"2x8x13", 165}, {"2x8x14", 176}, {"2x8x15", 188}, {"2x8x16", 200},
        {"2x9x9", 126}, {"2x9x10", 144}, {"2x9x11", 158}, {"2x9x12", 171}, {"2x9x13", 185}, {"2x9x14", 198},
        {"2x10x10", 155}, {"2x10x11", 174}, {"2x10x12", 188},
        {"2x11x11", 187},
        {"3x3x3", 23}, {"3x3x4", 29}, {"3x3x5", 36}, {"3x3x6", 42}, {"3x3x7", 49}, {"3x3x8", 56}, {"3x3x9", 63}, {"3x3x10", 71}, {"3x3x11", 78}, {"3x3x12", 84}, {"3x3x13", 91}, {"3x3x14", 98}, {"3x3x15", 105}, {"3x3x16", 112},
        {"3x4x4", 38}, {"3x4x5", 47}, {"3x4x6", 54}, {"3x4x7", 64}, {"3x4x8", 74}, {"3x4x9", 83}, {"3x4x10", 92}, {"3x4x11", 101}, {"3x4x12", 108}, {"3x4x13", 118}, {"3x4x14", 128}, {"3x4x15", 137}, {"3x4x16", 146},
        {"3x5x5", 58}, {"3x5x6", 70}, {"3x5x7", 81}, {"3x5x8", 92}, {"3x5x9", 105}, {"3x5x10", 115}, {"3x5x11", 128}, {"3x5x12", 139}, {"3x5x13", 150}, {"3x5x14", 162}, {"3x5x15", 173}, {"3x5x16", 184},
        {"3x6x6", 83}, {"3x6x7", 96}, {"3x6x8", 108}, {"3x6x9", 124}, {"3x6x10", 137}, {"3x6x11", 150}, {"3x6x12", 162}, {"3x6x13", 178}, {"3x6x14", 191}, {"3x6x15", 204}, {"3x6x16", 216},
        {"3x7x7", 113}, {"3x7x8", 128}, {"3x7x9", 145}, {"3x7x10", 160}, {"3x7x11", 177}, {"3x7x12", 192}, {"3x7x13", 209}, {"3x7x14", 224}, {"3x7x15", 241}, {"3x7x16", 256},
        {"3x8x8", 148}, {"3x8x9", 164}, {"3x8x10", 182}, {"3x8x11", 200}, {"3x8x12", 216}, {"3x8x13", 236}, {"3x8x14", 256}, {"3x8x15", 272}, {"3x8x16", 290},
        {"3x9x9", 187}, {"3x9x10", 207}, {"3x9x11", 227}, {"3x9x12", 246}, {"3x9x13", 268}, {"3x9x14", 288},
        {"3x10x10", 229}, {"3x10x11", 251}, {"3x10x12", 270},
        {"3x11x11", 278},
        {"4x4x4", 49}, {"4x4x5", 61}, {"4x4x6", 73}, {"4x4x7", 85}, {"4x4x8", 96}, {"4x4x9", 107}, {"4x4x10", 115}, {"4x4x11", 130}, {"4x4x12", 141}, {"4x4x13", 153}, {"4x4x14", 164}, {"4x4x15", 176}, {"4x4x16", 188},
        {"4x5x5", 76}, {"4x5x6", 90}, {"4x5x7", 104}, {"4x5x8", 118}, {"4x5x9", 132}, {"4x5x10", 146}, {"4x5x11", 160}, {"4x5x12", 175}, {"4x5x13", 192}, {"4x5x14", 207}, {"4x5x15", 221}, {"4x5x16", 236},
        {"4x6x6", 105}, {"4x6x7", 123}, {"4x6x8", 140}, {"4x6x9", 159}, {"4x6x10", 175}, {"4x6x11", 194}, {"4x6x12", 210}, {"4x6x13", 228}, {"4x6x14", 245}, {"4x6x15", 263}, {"4x6x16", 280},
        {"4x7x7", 144}, {"4x7x8", 164}, {"4x7x9", 187}, {"4x7x10", 207}, {"4x7x11", 227}, {"4x7x12", 246}, {"4x7x13", 267}, {"4x7x14", 285}, {"4x7x15", 307}, {"4x7x16", 324},
        {"4x8x8", 182}, {"4x8x9", 209}, {"4x8x10", 230}, {"4x8x11", 255}, {"4x8x12", 272}, {"4x8x13", 297}, {"4x8x14", 315}, {"4x8x15", 339}, {"4x8x16", 357},
        {"4x9x9", 225}, {"4x9x10", 255}, {"4x9x11", 279}, {"4x9x12", 300}, {"4x9x13", 332}, {"4x9x14", 357},
        {"4x10x10", 280}, {"4x10x11", 308}, {"4x10x12", 329},
        {"4x11x11", 342},
        {"5x5x5", 93}, {"5x5x6", 110}, {"5x5x7", 127}, {"5x5x8", 144}, {"5x5x9", 163}, {"5x5x10", 184}, {"5x5x11", 202}, {"5x5x12", 220}, {"5x5x13", 237}, {"5x5x14", 254}, {"5x5x15", 271}, {"5x5x16", 288},
        {"5x6x6", 130}, {"5x6x7", 150}, {"5x6x8", 170}, {"5x6x9", 197}, {"5x6x10", 217}, {"5x6x11", 240}, {"5x6x12", 258}, {"5x6x13", 280}, {"5x6x14", 300}, {"5x6x15", 320}, {"5x6x16", 340},
        {"5x7x7", 176}, {"5x7x8", 204}, {"5x7x9", 231}, {"5x7x10", 254}, {"5x7x11", 277}, {"5x7x12", 300}, {"5x7x13", 326}, {"5x7x14", 351}, {"5x7x15", 379}, {"5x7x16", 404},
        {"5x8x8", 230}, {"5x8x9", 262}, {"5x8x10", 287}, {"5x8x11", 313}, {"5x8x12", 333}, {"5x8x13", 365}, {"5x8x14", 391}, {"5x8x15", 423}, {"5x8x16", 451},
        {"5x9x9", 295}, {"5x9x10", 323}, {"5x9x11", 355}, {"5x9x12", 381}, {"5x9x13", 418}, {"5x9x14", 449},
        {"5x10x10", 352}, {"5x10x11", 390}, {"5x10x12", 421},
        {"5x11x11", 432},
        {"6x6x6", 153}, {"6x6x7", 183}, {"6x6x8", 203}, {"6x6x9", 225}, {"6x6x10", 252}, {"6x6x11", 276}, {"6x6x12", 294}, {"6x6x13", 322}, {"6x6x14", 343}, {"6x6x15", 371}, {"6x6x16", 392},
        {"6x7x7", 212}, {"6x7x8", 238}, {"6x7x9", 268}, {"6x7x10", 296}, {"6x7x11", 322}, {"6x7x12", 342}, {"6x7x13", 376}, {"6x7x14", 403}, {"6x7x15", 437}, {"6x7x16", 465},
        {"6x8x8", 266}, {"6x8x9", 296}, {"6x8x10", 329}, {"6x8x11", 357}, {"6x8x12", 378}, {"6x8x13", 418}, {"6x8x14", 448}, {"6x8x15", 486}, {"6x8x16", 518},
        {"6x9x9", 342}, {"6x9x10", 373}, {"6x9x11", 411}, {"6x9x12", 435}, {"6x9x13", 484}, {"6x9x14", 516},
        {"6x10x10", 406}, {"6x10x11", 454}, {"6x10x12", 490},
        {"6x11x11", 504},
        {"7x7x7", 250}, {"7x7x8", 279}, {"7x7x9", 316}, {"7x7x10", 346}, {"7x7x11", 378}, {"7x7x12", 404}, {"7x7x13", 443}, {"7x7x14", 475}, {"7x7x15", 513}, {"7x7x16", 548},
        {"7x8x8", 310}, {"7x8x9", 352}, {"7x8x10", 385}, {"7x8x11", 423}, {"7x8x12", 454}, {"7x8x13", 498}, {"7x8x14", 532}, {"7x8x15", 574}, {"7x8x16", 618},
        {"7x9x9", 399}, {"7x9x10", 437}, {"7x9x11", 482}, {"7x9x12", 520}, {"7x9x13", 567}, {"7x9x14", 604},
        {"7x10x10", 478}, {"7x10x11", 530}, {"7x10x12", 570},
        {"7x11x11", 584},
        {"8x8x8", 343}, {"8x8x9", 391}, {"8x8x10", 427}, {"8x8x11", 475}, {"8x8x12", 511}, {"8x8x13", 559}, {"8x8x14", 595}, {"8x8x15", 639}, {"8x8x16", 672},
        {"8x9x9", 435}, {"8x9x10", 487}, {"8x9x11", 539}, {"8x9x12", 570}, {"8x9x13", 631}, {"8x9x14", 671},
        {"8x10x10", 532}, {"8x10x11", 588}, {"8x10x12", 630},
        {"8x11x11", 646},
        {"9x9x9", 498}, {"9x9x10", 540}, {"9x9x11", 608}, {"9x9x12", 630}, {"9x9x13", 710}, {"9x9x14", 735},
        {"9x10x10", 600}, {"9x10x11", 662}, {"9x10x12", 705},
        {"9x11x11", 728},
        {"10x10x10", 651}, {"10x10x11", 719}, {"10x10x12", 770},
        {"10x11x11", 793},
        {"11x11x11", 873}
    };

    std::cout << "Initialized known ZT ranks" << std::endl;
}

template <typename Scheme>
void MetaFlipGraph<Scheme>::initializeKnownBinaryRanks() {
    dimension2knownRank = {
        {"2x2x2", 7}, {"2x2x3", 11}, {"2x2x4", 14}, {"2x2x5", 18}, {"2x2x6", 21}, {"2x2x7", 25}, {"2x2x8", 28}, {"2x2x9", 32}, {"2x2x10", 35}, {"2x2x11", 39}, {"2x2x12", 42}, {"2x2x13", 46}, {"2x2x14", 49}, {"2x2x15", 53}, {"2x2x16", 56},
        {"2x3x3", 15}, {"2x3x4", 20}, {"2x3x5", 25}, {"2x3x6", 30}, {"2x3x7", 35}, {"2x3x8", 40}, {"2x3x9", 45}, {"2x3x10", 50}, {"2x3x11", 55}, {"2x3x12", 60}, {"2x3x13", 65}, {"2x3x14", 70}, {"2x3x15", 75}, {"2x3x16", 80},
        {"2x4x4", 26}, {"2x4x5", 33}, {"2x4x6", 39}, {"2x4x7", 45}, {"2x4x8", 51}, {"2x4x9", 59}, {"2x4x10", 65}, {"2x4x11", 71}, {"2x4x12", 77}, {"2x4x13", 84}, {"2x4x14", 90}, {"2x4x15", 96}, {"2x4x16", 102},
        {"2x5x5", 40}, {"2x5x6", 47}, {"2x5x7", 55}, {"2x5x8", 63}, {"2x5x9", 72}, {"2x5x10", 80}, {"2x5x11", 87}, {"2x5x12", 94}, {"2x5x13", 102}, {"2x5x14", 110}, {"2x5x15", 118}, {"2x5x16", 127},
        {"2x6x6", 56}, {"2x6x7", 66}, {"2x6x8", 75}, {"2x6x9", 86}, {"2x6x10", 94}, {"2x6x11", 103}, {"2x6x12", 112}, {"2x6x13", 122}, {"2x6x14", 131}, {"2x6x15", 141}, {"2x6x16", 150},
        {"2x7x7", 76}, {"2x7x8", 88}, {"2x7x9", 100}, {"2x7x10", 110}, {"2x7x11", 121}, {"2x7x12", 131}, {"2x7x13", 142}, {"2x7x14", 152}, {"2x7x15", 164}, {"2x7x16", 175},
        {"2x8x8", 100}, {"2x8x9", 116}, {"2x8x10", 125}, {"2x8x11", 138}, {"2x8x12", 150}, {"2x8x13", 165}, {"2x8x14", 176}, {"2x8x15", 188}, {"2x8x16", 200},
        {"2x9x9", 126}, {"2x9x10", 140}, {"2x9x11", 154}, {"2x9x12", 168}, {"2x9x13", 185}, {"2x9x14", 198},
        {"2x10x10", 155}, {"2x10x11", 171}, {"2x10x12", 186},
        {"2x11x11", 187},
        {"3x3x3", 23}, {"3x3x4", 29}, {"3x3x5", 36}, {"3x3x6", 42}, {"3x3x7", 49}, {"3x3x8", 55}, {"3x3x9", 63}, {"3x3x10", 71}, {"3x3x11", 78}, {"3x3x12", 84}, {"3x3x13", 91}, {"3x3x14", 98}, {"3x3x15", 105}, {"3x3x16", 112},
        {"3x4x4", 38}, {"3x4x5", 47}, {"3x4x6", 54}, {"3x4x7", 64}, {"3x4x8", 73}, {"3x4x9", 83}, {"3x4x10", 92}, {"3x4x11", 101}, {"3x4x12", 108}, {"3x4x13", 118}, {"3x4x14", 127}, {"3x4x15", 137}, {"3x4x16", 146},
        {"3x5x5", 58}, {"3x5x6", 68}, {"3x5x7", 79}, {"3x5x8", 90}, {"3x5x9", 104}, {"3x5x10", 115}, {"3x5x11", 126}, {"3x5x12", 136}, {"3x5x13", 147}, {"3x5x14", 158}, {"3x5x15", 169}, {"3x5x16", 180},
        {"3x6x6", 83}, {"3x6x7", 96}, {"3x6x8", 108}, {"3x6x9", 122}, {"3x6x10", 136}, {"3x6x11", 150}, {"3x6x12", 162}, {"3x6x13", 178}, {"3x6x14", 191}, {"3x6x15", 204}, {"3x6x16", 216},
        {"3x7x7", 111}, {"3x7x8", 128}, {"3x7x9", 143}, {"3x7x10", 160}, {"3x7x11", 177}, {"3x7x12", 192}, {"3x7x13", 209}, {"3x7x14", 224}, {"3x7x15", 241}, {"3x7x16", 256},
        {"3x8x8", 145}, {"3x8x9", 164}, {"3x8x10", 180}, {"3x8x11", 200}, {"3x8x12", 216}, {"3x8x13", 236}, {"3x8x14", 256}, {"3x8x15", 270}, {"3x8x16", 290},
        {"3x9x9", 187}, {"3x9x10", 207}, {"3x9x11", 227}, {"3x9x12", 246}, {"3x9x13", 268}, {"3x9x14", 288},
        {"3x10x10", 229}, {"3x10x11", 251}, {"3x10x12", 270},
        {"3x11x11", 278},
        {"4x4x4", 47}, {"4x4x5", 60}, {"4x4x6", 73}, {"4x4x7", 85}, {"4x4x8", 94}, {"4x4x9", 107}, {"4x4x10", 115}, {"4x4x11", 130}, {"4x4x12", 141}, {"4x4x13", 153}, {"4x4x14", 164}, {"4x4x15", 176}, {"4x4x16", 188},
        {"4x5x5", 73}, {"4x5x6", 89}, {"4x5x7", 104}, {"4x5x8", 118}, {"4x5x9", 132}, {"4x5x10", 146}, {"4x5x11", 160}, {"4x5x12", 175}, {"4x5x13", 192}, {"4x5x14", 207}, {"4x5x15", 221}, {"4x5x16", 236},
        {"4x6x6", 105}, {"4x6x7", 123}, {"4x6x8", 140}, {"4x6x9", 159}, {"4x6x10", 175}, {"4x6x11", 194}, {"4x6x12", 210}, {"4x6x13", 228}, {"4x6x14", 245}, {"4x6x15", 263}, {"4x6x16", 280},
        {"4x7x7", 144}, {"4x7x8", 164}, {"4x7x9", 187}, {"4x7x10", 207}, {"4x7x11", 227}, {"4x7x12", 246}, {"4x7x13", 267}, {"4x7x14", 285}, {"4x7x15", 307}, {"4x7x16", 324},
        {"4x8x8", 182}, {"4x8x9", 209}, {"4x8x10", 230}, {"4x8x11", 255}, {"4x8x12", 272}, {"4x8x13", 297}, {"4x8x14", 315}, {"4x8x15", 339}, {"4x8x16", 357},
        {"4x9x9", 225}, {"4x9x10", 255}, {"4x9x11", 279}, {"4x9x12", 300}, {"4x9x13", 332}, {"4x9x14", 355},
        {"4x10x10", 280}, {"4x10x11", 308}, {"4x10x12", 329},
        {"4x11x11", 340},
        {"5x5x5", 93}, {"5x5x6", 110}, {"5x5x7", 127}, {"5x5x8", 144}, {"5x5x9", 163}, {"5x5x10", 183}, {"5x5x11", 200}, {"5x5x12", 217}, {"5x5x13", 237}, {"5x5x14", 254}, {"5x5x15", 271}, {"5x5x16", 288},
        {"5x6x6", 130}, {"5x6x7", 150}, {"5x6x8", 170}, {"5x6x9", 197}, {"5x6x10", 217}, {"5x6x11", 240}, {"5x6x12", 258}, {"5x6x13", 280}, {"5x6x14", 300}, {"5x6x15", 320}, {"5x6x16", 340},
        {"5x7x7", 176}, {"5x7x8", 204}, {"5x7x9", 229}, {"5x7x10", 254}, {"5x7x11", 277}, {"5x7x12", 300}, {"5x7x13", 326}, {"5x7x14", 351}, {"5x7x15", 379}, {"5x7x16", 404},
        {"5x8x8", 230}, {"5x8x9", 262}, {"5x8x10", 287}, {"5x8x11", 313}, {"5x8x12", 333}, {"5x8x13", 365}, {"5x8x14", 391}, {"5x8x15", 423}, {"5x8x16", 445},
        {"5x9x9", 295}, {"5x9x10", 323}, {"5x9x11", 355}, {"5x9x12", 381}, {"5x9x13", 418}, {"5x9x14", 449},
        {"5x10x10", 352}, {"5x10x11", 386}, {"5x10x12", 413},
        {"5x11x11", 432},
        {"6x6x6", 153}, {"6x6x7", 183}, {"6x6x8", 203}, {"6x6x9", 225}, {"6x6x10", 252}, {"6x6x11", 276}, {"6x6x12", 294}, {"6x6x13", 322}, {"6x6x14", 343}, {"6x6x15", 371}, {"6x6x16", 392},
        {"6x7x7", 212}, {"6x7x8", 238}, {"6x7x9", 268}, {"6x7x10", 296}, {"6x7x11", 322}, {"6x7x12", 342}, {"6x7x13", 376}, {"6x7x14", 403}, {"6x7x15", 437}, {"6x7x16", 465},
        {"6x8x8", 266}, {"6x8x9", 296}, {"6x8x10", 329}, {"6x8x11", 357}, {"6x8x12", 378}, {"6x8x13", 418}, {"6x8x14", 448}, {"6x8x15", 486}, {"6x8x16", 511},
        {"6x9x9", 342}, {"6x9x10", 373}, {"6x9x11", 411}, {"6x9x12", 435}, {"6x9x13", 484}, {"6x9x14", 516},
        {"6x10x10", 406}, {"6x10x11", 446}, {"6x10x12", 476},
        {"6x11x11", 504},
        {"7x7x7", 248}, {"7x7x8", 273}, {"7x7x9", 313}, {"7x7x10", 346}, {"7x7x11", 378}, {"7x7x12", 404}, {"7x7x13", 443}, {"7x7x14", 475}, {"7x7x15", 513}, {"7x7x16", 548},
        {"7x8x8", 302}, {"7x8x9", 352}, {"7x8x10", 385}, {"7x8x11", 423}, {"7x8x12", 454}, {"7x8x13", 498}, {"7x8x14", 532}, {"7x8x15", 574}, {"7x8x16", 618},
        {"7x9x9", 399}, {"7x9x10", 437}, {"7x9x11", 482}, {"7x9x12", 520}, {"7x9x13", 567}, {"7x9x14", 604},
        {"7x10x10", 478}, {"7x10x11", 526}, {"7x10x12", 564},
        {"7x11x11", 584},
        {"8x8x8", 329}, {"8x8x9", 391}, {"8x8x10", 427}, {"8x8x11", 475}, {"8x8x12", 511}, {"8x8x13", 559}, {"8x8x14", 595}, {"8x8x15", 639}, {"8x8x16", 672},
        {"8x9x9", 435}, {"8x9x10", 487}, {"8x9x11", 539}, {"8x9x12", 570}, {"8x9x13", 624}, {"8x9x14", 669},
        {"8x10x10", 532}, {"8x10x11", 588}, {"8x10x12", 630},
        {"8x11x11", 646},
        {"9x9x9", 498}, {"9x9x10", 540}, {"9x9x11", 608}, {"9x9x12", 630}, {"9x9x13", 710}, {"9x9x14", 735},
        {"9x10x10", 600}, {"9x10x11", 662}, {"9x10x12", 705},
        {"9x11x11", 728},
        {"10x10x10", 651}, {"10x10x11", 719}, {"10x10x12", 770},
        {"10x11x11", 793},
        {"11x11x11", 873}
    };

    std::cout << "Initialized known Z2 ranks" << std::endl;
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

    if (uniform(generator) < 0.5)
        scheme.swapSizes(generator);

    int index = generator() % schemesBest.size();
    bool resized = true;

    if (!scheme.tryMerge(schemesBest[index], generator, metaParameters.maxDimension, metaParameters.maxRank)) {
        double p = uniform(generator);

        if (p < 0.5) {
            resized = scheme.tryProject(generator, metaParameters.minDimension);
        }
        else {
            resized = scheme.tryExtend(generator, metaParameters.maxDimension, metaParameters.maxRank);
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
