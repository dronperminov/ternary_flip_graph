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

int BaseScheme::getDimension(int index) const {
    return dimension[index];
}

std::string BaseScheme::getDimension() const {
    std::stringstream ss;
    ss << dimension[0] << "x" << dimension[1] << "x" << dimension[2];
    return ss.str();
}

int BaseScheme::getAvailableFlips() const {
    return flips[0].size() + flips[1].size() + flips[2].size();
}

int BaseScheme::getAvailableFlips(int index) const {
    return flips[index].size();
}

double BaseScheme::getOmega() const {
    return 3 * log(rank) / log(dimension[0] * dimension[1] * dimension[2]);
}

FlipStructure BaseScheme::getOptimalStructure(std::mt19937 &generator, int iterations, double eps) const {
    FlipStructureOptimizer optimizer(dimension[0], dimension[1], dimension[2], rank);

    for (int i = 0; i < 3; i++)
        for (size_t j = 0; j < flips[i].size(); j++)
            optimizer.add(i, flips[i].index1(j), flips[i].index2(j));

    return optimizer.optimize(generator, iterations, eps);
}
