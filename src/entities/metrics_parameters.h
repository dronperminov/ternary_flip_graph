#pragma once

#include <iostream>
#include <string>
#include "arg_parser.h"

struct MetricsParameters {
    bool use;
    std::string path;

    void parse(const ArgParser &parser);
    friend std::ostream& operator<<(std::ostream& os, const MetricsParameters &metricsParameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
