#pragma once

#include <iostream>
#include "arg_parser.h"
#include "../utils.h"

struct FlipParameters {
    size_t flipIterations;
    size_t resetIterations;
    size_t minPlusIterations;
    size_t maxPlusIterations;
    int plusDiff;
    double sandwichingProbability;
    double reduceProbability;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const FlipParameters &flipParameters);
};
