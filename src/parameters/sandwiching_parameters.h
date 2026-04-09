#pragma once

#include <iostream>
#include <vector>
#include "../entities/arg_parser.h"
#include "../algebra/fraction.h"
#include "../utils.h"

struct SandwichingParameters {
    double probability;
    int maxDenominator;
    int minRows;
    int maxRows;
    int minColumns;
    int maxColumns;
    std::vector<Fraction> values;
    std::vector<Fraction> nonZero;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const SandwichingParameters &parameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
