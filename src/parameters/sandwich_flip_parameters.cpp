#include "sandwich_flip_parameters.h"

void SandwichFlipParameters::parse(const ArgParser &parser) {
    minSteps = std::stoi(parser["--min-steps"]);
    maxSteps = std::stoi(parser["--max-steps"]);

    minimizeNorm = parser.isSet("--minimize-norm");
    minimizeOmega = parser.isSet("--minimize-omega");
    maximizeFlips = parser.isSet("--maximize-flips");

    maxFractions = parser.isSet("--max-fractions") ? std::stoi(parser["--max-fractions"]) : INT_MAX;
    maxDenominator = parser.isSet("--max-denominator") ? std::stoi(parser["--max-denominator"]) : INT_MAX;
    maxNumerator = parser.isSet("--max-numerator") ? std::stoi(parser["--max-numerator"]) : INT_MAX;
    maxWeight = parser.isSet("--max-weight") ? std::stoi(parser["--max-weight"]) : INT_MAX;

    fixFractions = parser.isSet("--fix-fractions");
    check = parser["--weight-check"];
}

std::ostream& operator<<(std::ostream& os, const SandwichFlipParameters &parameters) {
    os << "SandwichFlip parameters:" << std::endl;
    os << "- steps: " << parameters.minSteps << " .. " << parameters.maxSteps << std::endl;
    os << "- minimize norm: " << (parameters.minimizeNorm ? "yes" : "no") << std::endl;
    os << "- minimize omega: " << (parameters.minimizeOmega ? "yes" : "no") << std::endl;
    os << "- maximize flips: " << (parameters.maximizeFlips ? "yes" : "no") << std::endl;
    os << "- fix fractions: " << (parameters.fixFractions ? "yes" : "no") << std::endl;
    os << "- check sequence: " << parameters.check << std::endl;

    if (parameters.maxFractions < INT_MAX || parameters.maxDenominator < INT_MAX || parameters.maxNumerator < INT_MAX || parameters.maxWeight < INT_MAX) {
        os << "- restrictions on parameters:" << std::endl;

        if (parameters.maxFractions < INT_MAX)
            os << "  - max fractions: " << parameters.maxFractions << std::endl;
        if (parameters.maxDenominator < INT_MAX)
            os << "  - max denominator: " << parameters.maxDenominator << std::endl;
        if (parameters.maxNumerator < INT_MAX)
            os << "  - max numerator: " << parameters.maxNumerator << std::endl;
        if (parameters.maxWeight < INT_MAX)
            os << "  - max weight: " << parameters.maxWeight << std::endl;
    }

    return os;
}

void SandwichFlipParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--min-steps", ArgType::Natural, "Minimum number of random steps", "1");
    parser.add("--max-steps", ArgType::Natural, "Maximum number of random steps", "5");

    parser.add("--maximize-flips", ArgType::Flag, "Check flips count during comparison first");
    parser.add("--minimize-omega", ArgType::Flag, "Check structure omega during comparison first");
    parser.add("--minimize-norm", ArgType::Flag, "Check flips count during comparison first");

    parser.add("--max-fractions", ArgType::Natural, "Maximum number of fractions");
    parser.add("--max-denominator", ArgType::Natural, "Maximum denominator of scheme");
    parser.add("--max-numerator", ArgType::Natural, "Maximum numerator of scheme");
    parser.add("--max-weight", ArgType::Natural, "Maximum weight of scheme");

    parser.add("--fix-fractions", ArgType::Flag, "Try to convert fractions to integers");
    parser.add("--weight-check", ArgType::String, "Sequence to check metrics", "dfnw");
}
