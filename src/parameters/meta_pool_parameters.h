#pragma once

#include <iostream>
#include "../entities/arg_parser.h"
#include "../utils.h"

struct MetaPoolParameters {
    bool use;
    size_t size;
    bool uniqueOnly;
    bool resume;
    bool liftOnly;
    double alternativesProbability;

    int mergeMaxDiff;
    int extendMaxDiff;
    int projectMaxDiff;
    int productMaxDiff;

    double mergeProbability;
    double extendProbability;
    double projectProbability;
    double productProbability;

    int projectMinN1;
    int projectMinN2;
    int projectMinN3;

    double selectRankScale;
    double metaRankScale;
    std::string prioritiesPath;

    void parse(const ArgParser &parser);
    void writeJSON(std::ostream &os) const;
    friend std::ostream& operator<<(std::ostream& os, const MetaPoolParameters &parameters);

    static void addToParser(ArgParser &parser, const std::string &sectionName);
};
