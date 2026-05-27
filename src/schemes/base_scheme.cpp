#include "base_scheme.h"

BaseScheme::BaseScheme() : boolDistribution(0, 1), ijkDistribution(0, 2) {
    for (int i = 0; i < 3; i++) {
        dimension[i] = 0;
        elements[i] = 0;
    }

    rank = 0;
}

int BaseScheme::getRank() const {
    return rank;
}

int BaseScheme::getCoefficientsCount() const {
    return rank * (elements[0] + elements[1] + elements[2]);
}

int BaseScheme::getDimension(int index) const {
    return dimension[index];
}

std::string BaseScheme::getDimension() const {
    std::stringstream ss;
    ss << dimension[0] << "x" << dimension[1] << "x" << dimension[2];
    return ss.str();
}

int BaseScheme::getElements(int index) const {
    return elements[index];
}

int BaseScheme::getAvailableFlips() const {
    return flips[0].size() + flips[1].size() + flips[2].size();
}

int BaseScheme::getAvailableFlips(int index) const {
    return flips[index].size();
}

int BaseScheme::getIndependentFlips() const {
    std::vector<int> used(rank, 0);

    for (int i = 0; i < 3; i++) {
        for (size_t j = 0; j < flips[i].size(); j++) {
            used[flips[i].index1(j)] |= 1 << i;
            used[flips[i].index2(j)] |= 1 << i;
        }
    }

    int independent = 0;

    for (int i = 0; i < rank; i++)
        if (used[i] && !(used[i] & (used[i] - 1)))
            independent++;

    return independent;
}

double BaseScheme::getOmega() const {
    return 3 * log(rank) / log(dimension[0] * dimension[1] * dimension[2]);
}

std::string BaseScheme::getStructureHash(std::mt19937 &generator) const {
    FlipStructureOptimizer optimizer = getStructureOptimizer();
    std::vector<std::vector<std::unordered_set<int>>> groups = optimizer.getGroups(generator);

    std::vector<std::vector<size_t>> groupSizes(3);
    for (int i = 0; i < 3; i++) {
        for (const auto &group : groups[i])
            groupSizes[i].push_back(group.size());

        std::sort(groupSizes[i].begin(), groupSizes[i].end());
    }

    std::stringstream ss;
    for (size_t size : groupSizes[0])
        ss << "u" << size;

    for (size_t size : groupSizes[1])
        ss << "v" << size;

    for (size_t size : groupSizes[2])
        ss << "w" << size;

    return ss.str();
}

FlipStructureOptimizer BaseScheme::getStructureOptimizer() const {
    FlipStructureOptimizer optimizer(dimension[0], dimension[1], dimension[2], rank);

    for (int i = 0; i < 3; i++)
        for (size_t j = 0; j < flips[i].size(); j++)
            optimizer.add(i, flips[i].index1(j), flips[i].index2(j));

    optimizer.preprocess();
    return optimizer;
}
