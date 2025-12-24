#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <random>

#include "../entities/flip_set.h"

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
};
