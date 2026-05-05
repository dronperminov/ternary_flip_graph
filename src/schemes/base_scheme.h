#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <cmath>

#include "../entities/flip_set.h"
#include "../entities/flip_structure_optimizer.h"

class BaseScheme {
protected:
    int dimension[3];
    int elements[3];
    int rank;

    FlipSet flips[3];

    std::uniform_int_distribution<int> boolDistribution;
    std::uniform_int_distribution<int> ijkDistribution;
public:
    BaseScheme();

    int getRank() const;
    int getCoefficientsCount() const;
    int getDimension(int index) const;
    std::string getDimension() const;
    int getAvailableFlips() const;
    int getAvailableFlips(int index) const;
    int getIndependentFlips() const;
    double getOmega() const;

    FlipStructureOptimizer getStructureOptimizer() const;
};
