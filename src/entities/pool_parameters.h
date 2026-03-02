#pragma once

#include <iostream>
#include "arg_parser.h"
#include "../utils.h"

struct PoolParameters {
    bool use;
    size_t maxIterations;
    size_t size;
    size_t minSize;
    std::string selectStrategy;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const PoolParameters &poolParameters);
};
