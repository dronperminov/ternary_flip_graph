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

void MetricsParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--save-metrics", ArgType::Flag, "Evaluate and save metrics");
    parser.add("--metrics-path", ArgType::Path, "Path to file with metrics", "schemes/metrics.jsonl");
}
