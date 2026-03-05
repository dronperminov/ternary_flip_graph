#include "pool_parameters.h"

void PoolParameters::parse(const ArgParser &parser) {
    use = parser.isSet("--use-pool");
    size = parseNatural(parser["--pool-size"]);
    minSize = parseNatural(parser["--pool-min-size"]);
    maxIterations = parseNatural(parser["--pool-max-iterations"]);
    selectStrategy = parser["--pool-select-strategy"];
}

std::ostream& operator<<(std::ostream& os, const PoolParameters &poolParameters) {
    if (poolParameters.use) {
        os << "Pool parameters:" << std::endl;
        os << "- size: " << poolParameters.size << " (minimal: " << poolParameters.minSize << ")" << std::endl;
        os << "- max iterations: " << poolParameters.maxIterations << std::endl;
        os << "- select strategy: " << poolParameters.selectStrategy << std::endl;
    }

    return os;
}

void PoolParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--use-pool", ArgType::Flag, "Use pool strategy");
    parser.add("--pool-size", ArgType::Natural, "Optimal size of pool", "1K");
    parser.add("--pool-min-size", ArgType::Natural, "Minimal size of pool", "5");
    parser.add("--pool-max-iterations", ArgType::Natural, "Max random walk iterations to reach min pool size", "1K");
    parser.addChoices("--pool-select-strategy", ArgType::String, "Pool selection strategy", {"uniform", "flips"}, "uniform");
}
