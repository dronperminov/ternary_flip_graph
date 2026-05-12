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

std::string BaseScheme::getStructureHash() const {
    std::vector<std::vector<int>> counts(rank, std::vector<int>(3, 0));
    for (int i = 0; i < 3; i++) {
        for (size_t j = 0; j < flips[i].size(); j++) {
            counts[flips[i].index1(j)][i]++;
            counts[flips[i].index2(j)][i]++;
        }
    }

    int max = 0;
    for (int i = 0; i < rank; i++)
        for (int j = 0; j < 3; j++)
            max = std::max(max, counts[i][j]);

    int digits = digitsCount(max);
    std::vector<std::string> lines;

    for (int i = 0; i < rank; i++) {
        std::stringstream line;
        line << std::setw(digits) << std::setfill('0') << counts[i][0] << "-";
        line << std::setw(digits) << std::setfill('0') << counts[i][1] << "-";
        line << std::setw(digits) << std::setfill('0') << counts[i][2];
        lines.push_back(line.str());
    }

    std::sort(lines.begin(), lines.end());

    std::stringstream hash;
    for (const std::string &line : lines)
        hash << line << "_";
    return hash.str();
}

FlipStructureOptimizer BaseScheme::getStructureOptimizer() const {
    FlipStructureOptimizer optimizer(dimension[0], dimension[1], dimension[2], rank);

    for (int i = 0; i < 3; i++)
        for (size_t j = 0; j < flips[i].size(); j++)
            optimizer.add(i, flips[i].index1(j), flips[i].index2(j));

    optimizer.preprocess();
    return optimizer;
}
