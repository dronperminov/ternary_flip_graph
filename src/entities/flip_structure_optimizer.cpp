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

void FlipStructureOptimizer::preprocess() {
    std::vector<int> used(rank, 0);
    for (const Flip &flip : flips) {
        used[flip.i] |= 1 << flip.p;
        used[flip.j] |= 1 << flip.p;
    }

    independentFlips.clear();
    dependentFlips.clear();

    for (const Flip &flip : flips) {
        if (isPowerOfTwo(used[flip.i]) && isPowerOfTwo(used[flip.j]))
            independentFlips.push_back(flip);
        else
            dependentFlips.push_back(flip);
    }
}

std::vector<Flip> FlipStructureOptimizer::selectRandomFlips(const std::vector<Flip> &flips, std::mt19937 &generator) const {
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
    std::vector<Flip> selected = selectRandomFlips(flips, generator);
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
    double n = dimension[0];
    double m = dimension[1];
    double p = dimension[2];

    double score = -std::pow(n * m * p, omega);

    for (const auto &a : structure) {
        for (const auto &b : structure) {
            for (const auto &c : structure) {
                double as = a.s;
                double bs = b.s;
                double cs = c.s;

                score += as * bs * cs * std::pow(a.n * b.m * c.p, omega - 2) * b.n * c.n * a.m * c.m * a.p * b.p;
            }
        }
    }

    return score;
}

double FlipStructureOptimizer::df(double omega, const std::vector<FlipStructureNode> &structure) const {
    double n = dimension[0];
    double m = dimension[1];
    double p = dimension[2];

    double score = -std::log(n * m * p) * std::pow(n * m * p, omega);

    for (const auto &a : structure) {
        for (const auto &b : structure) {
            for (const auto &c : structure) {
                double as = a.s;
                double bs = b.s;
                double cs = c.s;

                score += as * bs * cs * std::log(a.n * b.m * c.p) * std::pow(a.n * b.m * c.p, omega - 2) * b.n * c.n * a.m * c.m * a.p * b.p;
            }
        }
    }

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

std::vector<Flip> FlipStructureOptimizer::getFlips() const {
    return flips;
}

std::unordered_map<std::string, int> FlipStructureOptimizer::getSerendipitousRanks(std::mt19937 &generator, const std::unordered_map<std::string, int> &dimension2rank, int iterations, int maxN) const {
    std::unordered_map<std::string, int> dimension2serendipitousRank;

    for (int n1 = 1; n1 <= maxN; n1++)
        for (int n2 = n1; n2 <= maxN; n2++)
            for (int n3 = n2; n3 <= maxN; n3++)
                dimension2serendipitousRank[getDimension(n1, n2, n3)] = n1*n2*n3;

    if (dependentFlips.empty())
        iterations = 1;

    for (int iteration = 0; iteration < iterations; iteration++) {
        std::vector<Flip> selected = selectRandomFlips(dependentFlips, generator);
        for (const Flip &flip : independentFlips)
            selected.push_back(flip);

        std::vector<std::unordered_set<int>> u = groupFlips(selected, 0);
        std::vector<std::unordered_set<int>> v = groupFlips(selected, 1);
        std::vector<std::unordered_set<int>> w = groupFlips(selected, 2);
        std::unordered_set<int> flipIndices;

        for (const Flip& flip: selected) {
            flipIndices.insert(flip.i);
            flipIndices.insert(flip.j);
        }

        for (int n1 = 1; n1 <= maxN && dimension[0] * n1 <= maxN; n1++) {
            for (int n2 = 1; n2 <= maxN && dimension[1] * n2 <= maxN; n2++) {
                for (int n3 = 1; n3 <= maxN && dimension[2] * n3 <= maxN; n3++) {
                    if (n1 == 1 && n2 == 1 && n3 == 1)
                        continue;

                    int serendipitousRank = (rank - flipIndices.size()) * dimension2rank.at(getDimension(n1, n2, n3, true));

                    for (const auto& group : u) {
                        auto it = dimension2rank.find(getDimension(n1, n2, n3 * group.size(), true));
                        serendipitousRank += it == dimension2rank.end() ? n1 * n2 * n3 * group.size() : it->second;
                    }

                    for (const auto& group : v) {
                        auto it = dimension2rank.find(getDimension(n1 * group.size(), n2, n3, true));
                        serendipitousRank += it == dimension2rank.end() ? n1 * n2 * n3 * group.size() : it->second;
                    }

                    for (const auto& group : w) {
                        auto it = dimension2rank.find(getDimension(n1, n2 * group.size(), n3, true));
                        serendipitousRank += it == dimension2rank.end() ? n1 * n2 * n3 * group.size() : it->second;
                    }

                    std::string key = getDimension(dimension[0] * n1, dimension[1] * n2, dimension[2] * n3, true);
                    dimension2serendipitousRank[key] = std::min(dimension2serendipitousRank.at(key), serendipitousRank);
                }
            }
        }
    }

    return dimension2serendipitousRank;
}

std::vector<std::vector<std::unordered_set<int>>> FlipStructureOptimizer::getGroups(std::mt19937 &generator) const {
    std::vector<Flip> selected = selectRandomFlips(dependentFlips, generator);
    for (const Flip &flip : independentFlips)
        selected.push_back(flip);

    std::vector<std::vector<std::unordered_set<int>>> groups;
    for (int i = 0; i < 3; i++)
        groups.push_back(groupFlips(selected, i));

    return groups;
}

std::ostream& operator<<(std::ostream &os, const std::vector<FlipStructureNode> &structure) {
    for (size_t i = 0; i < structure.size(); i++) {
        FlipStructureNode node = structure[i];
        os << (i > 0 ? "+" : "") << node.s << "<" << node.n << "," << node.m << "," << node.p << ">";
    }

    return os;
}

std::ostream& operator<<(std::ostream &os, const FlipStructure &structure) {
    os << "[";

    for (size_t i = 0; i < structure.structure.size(); i++) {
        FlipStructureNode node = structure.structure[i];
        os << (i > 0 ? ", " : "") << "[" << node.s << ", [" << node.n << ", " << node.m << ", " << node.p << "]]";
    }

    os << "]";
    return os;
}
