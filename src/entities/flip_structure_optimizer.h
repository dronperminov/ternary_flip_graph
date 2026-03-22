#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

struct Flip {
    int p;
    int i;
    int j;
};

struct FlipStructureNode {
    int s;
    int n;
    int m;
    int p;
};

struct FlipStructure {
    double omega;
    std::vector<FlipStructureNode> structure;
};

class FlipStructureOptimizer {
    int dimension[3];
    int rank;
    std::vector<Flip> flips;

    std::vector<Flip> selectRandomFlips(std::mt19937 &generator) const;
    std::vector<std::unordered_set<int>> groupFlips(const std::vector<Flip> &flips, int p) const;
    std::vector<int> countSizes(const std::vector<std::unordered_set<int>> &components) const;
    std::vector<FlipStructureNode> selectRandomStructure(std::mt19937 &generator) const;

    double f(double omega, const std::vector<FlipStructureNode> &structure) const;
    double df(double omega, const std::vector<FlipStructureNode> &structure) const;
    double findOmega(const std::vector<FlipStructureNode> &structure, double a, double b, double eps) const;
    double findOmega(const std::vector<FlipStructureNode> &structure, double x0, double eps) const;
public:
    FlipStructureOptimizer(int n1, int n2, int n3, int rank);

    void add(int p, int i, int j);

    FlipStructure optimize(std::mt19937 &generator, int iterations, double eps);
};
