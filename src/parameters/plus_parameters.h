#pragma once

#include <iostream>
#include <vector>
#include "../entities/arg_parser.h"
#include "../algebra/fraction.h"
#include "../utils.h"

struct PlusParameters {
    double probability;
    int iterations;
    std::vector<Fraction> values;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const PlusParameters &parameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
