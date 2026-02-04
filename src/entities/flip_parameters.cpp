#include "flip_parameters.h"

void FlipParameters::parse(const ArgParser &parser) {
    flipIterations = parseNatural(parser["--flip-iterations"]);
    minPlusIterations = parseNatural(parser["--min-plus-iterations"]);
    maxPlusIterations = parseNatural(parser["--max-plus-iterations"]);
    resetIterations = parseNatural(parser["--reset-iterations"]);
    plusDiff = std::stoi(parser["--plus-diff"]);
    sandwichingProbability = std::stod(parser["--sandwiching-probability"]);
    reduceProbability = std::stod(parser["--reduce-probability"]);
}

std::ostream& operator<<(std::ostream& os, const FlipParameters &flipParameters) {
    os << "Random walk parameters:" << std::endl;
    os << "- flip iterations: " << prettyInt(flipParameters.flipIterations) << std::endl;
    os << "- plus iterations: " << prettyInt(flipParameters.minPlusIterations) << " .. " << prettyInt(flipParameters.maxPlusIterations) << std::endl;
    os << "- reset iterations: " << prettyInt(flipParameters.resetIterations) << std::endl;
    os << "- plus diff: " << flipParameters.plusDiff << std::endl;
    os << "- sandwiching probability: " << flipParameters.sandwichingProbability << std::endl;
    os << "- reduce probability: " << flipParameters.reduceProbability << std::endl;
    return os;
}
