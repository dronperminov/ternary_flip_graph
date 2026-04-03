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
    int getDimension(int index) const;
    std::string getDimension() const;
    int getAvailableFlips() const;
    int getAvailableFlips(int index) const;
    double getOmega() const;

    FlipStructure getOptimalStructure(std::mt19937 &generator, int iterations = 250, double eps = 1e-15) const;
};
