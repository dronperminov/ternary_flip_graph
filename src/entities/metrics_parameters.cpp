#include "metrics_parameters.h"

void MetricsParameters::parse(const ArgParser &parser) {
    use = parser.isSet("--save-metrics");
    path = parser["--metrics-path"];
}

std::ostream& operator<<(std::ostream& os, const MetricsParameters &metricsParameters) {
    if (metricsParameters.use) {
        os << "Metrics parameters:" << std::endl;
        os << "- path: " << metricsParameters.path << std::endl;
    }

    return os;
}
