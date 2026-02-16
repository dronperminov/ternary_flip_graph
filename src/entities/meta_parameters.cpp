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
