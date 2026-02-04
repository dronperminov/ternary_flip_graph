#include "meta_parameters.h"

void MetaParameters::parse(const ArgParser &parser) {
    probability = std::stod(parser["--meta-probability"]);
    minDimension = std::stoi(parser["--meta-min-dimension"]);
    maxDimension = std::stoi(parser["--meta-max-dimension"]);
    maxRank = std::stoi(parser["--meta-max-rank"]);
}

std::ostream& operator<<(std::ostream& os, const MetaParameters &metaParameters) {
    os << "Meta operations parameters:" << std::endl;
    os << "- meta probability: " << metaParameters.probability << std::endl;
    os << "- meta dimensions: " << metaParameters.minDimension << " .. " << metaParameters.maxDimension << std::endl;
    os << "- meta max rank: " << metaParameters.maxRank << std::endl;
    return os;
}
