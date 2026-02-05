#pragma once

#include <iostream>
#include "arg_parser.h"
#include "../utils.h"

struct MetaParameters {
    double probability;
    std::string strategy;
    int minDimension;
    int maxDimension;
    int maxRank;
    int maxRankDiff;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const MetaParameters &metaParameters);
};
