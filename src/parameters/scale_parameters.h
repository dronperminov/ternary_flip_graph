#pragma once

#include <iostream>
#include <vector>
#include "../entities/arg_parser.h"
#include "../algebra/fraction.h"
#include "../utils.h"

struct ScaleParameters {
    double probability;
    double fullProbability;
    std::vector<Fraction> values;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const ScaleParameters &parameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
