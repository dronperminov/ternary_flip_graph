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
    int flipIterations;
    int plusIterations;
    double reduceProbability;
    int seed;
    int topCount;

    std::vector<TernaryScheme> schemes;
    std::vector<TernaryScheme> schemesBest;
    std::vector<int> flips;
    std::vector<int> bestRanks;
    std::vector<int> indices;
    int bestRank;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
public:
    FlipGraph(int count, const std::string outputPath, int threads, int flipIterations, int plusIterations, double reduceProbability, int seed, int topCount);

    void run(const TernaryScheme &scheme, int targetRank);
private:
    void initialize(const TernaryScheme &scheme);
    void runIteration();
    void updateBest(int iteration);
    void report(int iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;
    void randomWalk(TernaryScheme &scheme, TernaryScheme &schemeBest, int &flipsCount, int &bestRank, std::mt19937 &generator);

    bool compare(int index1, int index2) const;
    std::string getSavePath(const TernaryScheme &scheme, int iteration, const std::string path) const;
    std::string prettyInt(int value) const;
    std::string prettyDimension(const TernaryScheme &scheme) const;
    std::string prettyTime(double elapsed) const;
};
