#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <omp.h>

#include "../schemes/ternary_scheme.h"


class FlipGraph {
    int count;
    std::string outputPath;
    int threads;
    size_t flipIterations;
    size_t plusIterations;
    size_t resetIterations;
    double reduceProbability;
    int seed;
    int topCount;

    std::vector<TernaryScheme> schemes;
    std::vector<TernaryScheme> schemesBest;
    std::vector<size_t> flips;
    std::vector<int> bestRanks;
    std::vector<int> indices;
    int bestRank;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
public:
    FlipGraph(int count, const std::string outputPath, int threads, size_t flipIterations, size_t plusIterations, double reduceProbability, int seed, int topCount);

    void run(const TernaryScheme &scheme, int targetRank);
private:
    void initialize(const TernaryScheme &scheme);
    void runIteration();
    void updateBest(size_t iteration);
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void randomWalk(TernaryScheme &scheme, TernaryScheme &schemeBest, size_t &flipsCount, int &bestRank, std::mt19937 &generator);

    bool compare(int index1, int index2) const;
    std::string getSavePath(const TernaryScheme &scheme, int iteration, const std::string path) const;
    std::string prettyInt(size_t value) const;
    std::string prettyDimension(const TernaryScheme &scheme) const;
    std::string prettyTime(double elapsed) const;
};
