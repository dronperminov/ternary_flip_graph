#pragma once

#include <iostream>
#include "../entities/arg_parser.h"
#include "../utils.h"

struct PoolParameters {
    bool use;
    size_t maxIterations;
    size_t size;
    size_t minSize;
    bool uniqueOnly;
    std::string selectStrategy;

    void parse(const ArgParser &parser);
    void writeJSON(std::ostream &os) const;
    friend std::ostream& operator<<(std::ostream& os, const PoolParameters &poolParameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
