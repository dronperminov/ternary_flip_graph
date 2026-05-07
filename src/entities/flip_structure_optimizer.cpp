#include "flip_structure_optimizer.h"

FlipStructureOptimizer::FlipStructureOptimizer(int n1, int n2, int n3, int rank) {
    dimension[0] = n1;
    dimension[1] = n2;
    dimension[2] = n3;
    this->rank = rank;
}

void FlipStructureOptimizer::add(int p, int i, int j) {
    flips.push_back({p, i, j});
}

std::vector<Flip> FlipStructureOptimizer::selectRandomFlips(std::mt19937 &generator) const {
    std::vector<Flip> available(flips);
    std::unordered_set<int> ignored[3];
    std::vector<Flip> selected;

    while (!available.empty()) {
        Flip flip = available[generator() % available.size()];

        for (int p = 0; p < 3; p++) {
            if (p != flip.p) {
                ignored[p].insert(flip.i);
                ignored[p].insert(flip.j);
            }
        }

        auto it = std::remove_if(available.begin(), available.end(), [ignored, flip](const Flip& f) {
            return ignored[f.p].find(f.i) != ignored[f.p].end() || ignored[f.p].find(f.j) != ignored[f.p].end() || (f.p == flip.p && f.i == flip.i && f.j == flip.j);
        });

        available.erase(it, available.end());
        selected.push_back(flip);
    }

    return selected;
}

std::vector<std::unordered_set<int>> FlipStructureOptimizer::groupFlips(const std::vector<Flip> &flips, int p) const {
    std::unordered_map<int, int> index2component;
    std::vector<std::unordered_set<int>> components;

    for (const auto &flip : flips) {
        if (flip.p != p)
            continue;

        if (index2component.find(flip.i) == index2component.end()) {
            index2component[flip.i] = components.size();
            components.push_back({flip.i});
        }

        if (index2component.find(flip.j) == index2component.end()) {
            index2component[flip.j] = components.size();
            components.push_back({flip.j});
        }
    }

    for (const auto &flip : flips) {
        if (flip.p != p)
            continue;

        int ci = index2component[flip.i];
        int cj = index2component[flip.j];
        if (ci == cj)
            continue;

        for (int index : components[cj]) {
            index2component[index] = ci;
            components[ci].insert(index);
        }

        components[cj].clear();
    }

    auto it = std::remove_if(components.begin(), components.end(), [](const std::unordered_set<int> &component) {
        return component.empty();
    });
    components.erase(it, components.end());
    return components;
}

std::vector<int> FlipStructureOptimizer::countSizes(const std::vector<std::unordered_set<int>> &components) const {
    size_t maxCount = 2;

    for (const auto &component : components)
        maxCount = std::max(maxCount, component.size());

    std::vector<int> counts(maxCount + 1, 0);

    for (const auto &component : components)
        counts[component.size()]++;

    return counts;
}

std::vector<FlipStructureNode> FlipStructureOptimizer::selectRandomStructure(std::mt19937 &generator) const {
    std::vector<Flip> selected = selectRandomFlips(generator);
    std::vector<std::unordered_set<int>> u = groupFlips(selected, 0);
    std::vector<std::unordered_set<int>> v = groupFlips(selected, 1);
    std::vector<std::unordered_set<int>> w = groupFlips(selected, 2);

    std::vector<int> sizesU = countSizes(u);
    std::vector<int> sizesV = countSizes(v);
    std::vector<int> sizesW = countSizes(w);

    std::vector<FlipStructureNode> structure;

    for (size_t size = 0; size < sizesU.size(); size++)
        if (sizesU[size])
            structure.push_back({sizesU[size], 1, 1, (int)size});

    for (size_t size = 0; size < sizesV.size(); size++)
        if (sizesV[size])
            structure.push_back({sizesV[size], (int)size, 1, 1});

    for (size_t size = 0; size < sizesW.size(); size++)
        if (sizesW[size])
            structure.push_back({sizesW[size], 1, (int)size, 1});

    int total = 0;
    for (const auto &node : structure)
        total += node.s * node.n * node.m * node.p;

    if (total < rank)
        structure.push_back({rank - total, 1, 1, 1});

    return structure;
}

double FlipStructureOptimizer::f(double omega, const std::vector<FlipStructureNode> &structure) const {
    int n = dimension[0];
    int m = dimension[1];
    int p = dimension[2];

    double score = -std::pow(n * m * p, omega);

    for (const auto &a : structure)
        for (const auto &b : structure)
            for (const auto &c : structure)
                score += a.s * b.s * c.s * std::pow(a.n * b.m * c.p, omega - 2) * b.n * c.n * a.m * c.m * a.p * b.p;

    return score;
}

double FlipStructureOptimizer::df(double omega, const std::vector<FlipStructureNode> &structure) const {
    int n = dimension[0];
    int m = dimension[1];
    int p = dimension[2];

    double score = -std::log(n * m * p) * std::pow(n * m * p, omega);

    for (const auto &a : structure)
        for (const auto &b : structure)
            for (const auto &c : structure)
                score += a.s * b.s * c.s * std::log(a.n * b.m * c.p) * std::pow(a.n * b.m * c.p, omega - 2) * b.n * c.n * a.m * c.m * a.p * b.p;

    return score;
}

double FlipStructureOptimizer::findOmega(const std::vector<FlipStructureNode> &structure, double a, double b, double eps) const {
    double fa = f(a, structure);

    while (b - a > eps) {
        double x = (a + b) / 2;
        double fx = f(x, structure);

        if (fa * fx < 0) {
            b = x;
        }
        else {
            a = x;
            fa = fx;
        }
    }

    return (a + b) / 2;
}

double FlipStructureOptimizer::findOmega(const std::vector<FlipStructureNode> &structure, double x0, double eps) const {
    double x1 = x0 - f(x0, structure) / df(x0, structure);

    while (fabs(x1 - x0) > eps) {
        x0 = x1;
        x1 = x0 - f(x0, structure) / df(x0, structure);
    }

    return x1;
}

FlipStructure FlipStructureOptimizer::optimize(std::mt19937 &generator, int iterations, double eps) {
    FlipStructure optimized;
    optimized.omega = 3;
    optimized.structure = {};

    for (int iteration = 0; iteration < iterations; iteration++) {
        std::vector<FlipStructureNode> structure = selectRandomStructure(generator);

        double omega = findOmega(structure, 2.8, eps);
        if (omega < optimized.omega) {
            optimized.omega = omega;
            optimized.structure = structure;
        }
    }

    return optimized;
}

int FlipStructureOptimizer::getSerendipitousRank(std::mt19937 &generator, int dimension[3], const std::unordered_map<std::string, int> &dimension2rank, int iterations) {
    std::string dim = getDimension(dimension[0], dimension[1], dimension[2], true);
    if (dimension2rank.find(dim) == dimension2rank.end())
        return -1;

    int bestRank = 10000000;

    for (int iteration = 0; iteration < iterations; iteration++) {
        std::vector<Flip> selected = selectRandomFlips(generator);
        std::vector<std::unordered_set<int>> u = groupFlips(selected, 0);
        std::vector<std::unordered_set<int>> v = groupFlips(selected, 1);
        std::vector<std::unordered_set<int>> w = groupFlips(selected, 2);
        std::unordered_set<int> flipIndices;

        for (const Flip& flip: selected) {
            flipIndices.insert(flip.i);
            flipIndices.insert(flip.j);
        }

        int serendipitousRank = (rank - flipIndices.size()) * dimension2rank.at(dim);

        for (const auto& group : u) {
            std::string groupDimension = getDimension(dimension[0], dimension[1], dimension[2] * group.size(), true);
            if (dimension2rank.find(groupDimension) == dimension2rank.end())
                return -1;

            serendipitousRank += dimension2rank.at(groupDimension);
        }

        for (const auto& group : v) {
            std::string groupDimension = getDimension(dimension[0] * group.size(), dimension[1], dimension[2], true);
            if (dimension2rank.find(groupDimension) == dimension2rank.end())
                return -1;

            serendipitousRank += dimension2rank.at(groupDimension);
        }

        for (const auto& group : w) {
            std::string groupDimension = getDimension(dimension[0], dimension[1] * group.size(), dimension[2], true);
            if (dimension2rank.find(groupDimension) == dimension2rank.end())
                return -1;

            serendipitousRank += dimension2rank.at(groupDimension);
        }

        if (serendipitousRank < bestRank)
            bestRank = serendipitousRank;
    }

    return bestRank;
}
