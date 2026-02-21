#pragma once

#include <iostream>
#include "arg_parser.h"
#include "../utils.h"

struct PoolParameters {
    bool use;
    size_t size;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const PoolParameters &poolParameters);
};
