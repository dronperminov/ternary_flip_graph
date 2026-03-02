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
