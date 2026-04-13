#include "meta_parameters.h"

void MetaParameters::parse(const ArgParser &parser) {
    probability = std::stod(parser["--meta-probability"]);
    strategy = parser["--meta-strategy"];
    minDimension = std::stoi(parser["--meta-min-dimension"]);
    maxDimension = std::stoi(parser["--meta-max-dimension"]);
    maxRank = std::stoi(parser["--meta-max-rank"]);
    maxRankDiff = std::stoi(parser["--meta-max-rank-diff"]);
}

std::ostream& operator<<(std::ostream& os, const MetaParameters &metaParameters) {
    if (metaParameters.probability == 0)
        return os << "Meta operations: not used" << std::endl;

    os << "Meta operations parameters:" << std::endl;
    os << "- meta probability: " << metaParameters.probability << std::endl;
    os << "- meta strategy: " << metaParameters.strategy << std::endl;
    os << "- meta dimensions: " << metaParameters.minDimension << " .. " << metaParameters.maxDimension << std::endl;
    os << "- meta max rank: " << metaParameters.maxRank << std::endl;
    os << "- meta max rank diff: " << metaParameters.maxRankDiff << std::endl;
    return os;
}

void MetaParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--meta-probability", ArgType::Real, "Probability of call meta operations, from 0.0 to 1.0", "0");
    parser.addChoices("--meta-strategy", ArgType::String, "Strategy of meta operations", {"default", "proj", "ext"}, "default");
    parser.add("--meta-min-dimension", ArgType::Natural, "Min dimension for project meta operation", "2");
    parser.add("--meta-max-dimension", ArgType::Natural, "Max dimension for merge/extend meta operations", "16");
    parser.add("--meta-max-rank", ArgType::Natural, "Max rank for merge/extend meta operations", "350");
    parser.add("--meta-max-rank-diff", ArgType::UInt, "Max rank difference for reset to initial", "10");
}
