#pragma once

#include <iostream>
#include <string>
#include <climits>
#include "../entities/arg_parser.h"
#include "../utils.h"

struct SandwichFlipParameters {
    int minSteps;
    int maxSteps;
    bool minimizeNorm;
    bool minimizeOmega;
    bool maximizeFlips;

    int maxFractions;
    int maxDenominator;
    int maxNumerator;
    int maxWeight;

    bool fixFractions;
    std::string check;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const SandwichFlipParameters &parameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
