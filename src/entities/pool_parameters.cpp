#include "pool_parameters.h"

void PoolParameters::parse(const ArgParser &parser) {
    use = parser.isSet("--use-pool");
    size = parseNatural(parser["--pool-size"]);
}

std::ostream& operator<<(std::ostream& os, const PoolParameters &poolParameters) {
    if (poolParameters.use) {
        os << "Pool parameters:" << std::endl;
        os << "- size: " << poolParameters.size << std::endl;
    }

    return os;
}
