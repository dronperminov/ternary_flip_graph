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

void FlipParameters::writeJSON(std::ostream &os) const {
    os << "{";
    os << "\"iterations\": " << flipIterations << ", ";
    os << "\"min_plus_iterations\": " << minPlusIterations << ", ";
    os << "\"max_plus_iterations\": " << maxPlusIterations << ", ";
    os << "\"reset_iterations\": " << resetIterations << ", ";
    os << "\"plus_diff\": " << plusDiff << ", ";
    os << "\"sandwiching_probability\": " << sandwichingProbability << ", ";
    os << "\"reduce_probability\": " << reduceProbability;
    os << "}";
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

void FlipParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--flip-iterations", ArgType::Natural, "Flip iterations before reporting", "1M");
    parser.add("--min-plus-iterations", ArgType::Natural, "Minimum period for plus operator calls", "5K");
    parser.add("--max-plus-iterations", ArgType::Natural, "Maximum period for plus operator calls", "100K");
    parser.add("--reset-iterations", ArgType::Natural, "Total iterations before reset", "10B");
    parser.add("--plus-diff", ArgType::UInt, "Maximum rank difference for plus operations", "4");
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0");
    parser.add("--reduce-probability", ArgType::Real, "Probability of reduce operation, from 0.0 to 1.0", "0");
}
