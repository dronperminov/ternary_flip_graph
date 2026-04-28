#include "scale_parameters.h"

void ScaleParameters::parse(const ArgParser &parser) {
    maxRows = std::stoi(parser["--scale-max-rows"]);
    probability = std::stod(parser["--scale-probability"]);
    fullProbability = std::stod(parser["--scale-full-probability"]);
    values = parseFractions(parser["--scale-values"]);
}

std::ostream& operator<<(std::ostream& os, const ScaleParameters &parameters) {
    if (parameters.probability == 0) {
        os << "Scale operation: not used" << std::endl;
        return os;
    }

    os << "Scale parameters:" << std::endl;
    os << "- max rows: " << parameters.maxRows << std::endl;
    os << "- probability: " << parameters.probability << std::endl;
    os << "- full probability: " << parameters.fullProbability << std::endl;
    os << "- values: [";
    for (size_t i = 0; i < parameters.values.size(); i++)
        os << (i > 0 ? ", " : "") << parameters.values[i].pretty();
    os << "]" << std::endl;
    return os;
}

void ScaleParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--scale-max-rows", ArgType::Natural, "Scale max rows", "1");
    parser.add("--scale-probability", ArgType::Real, "Probability of scale operation, from 0.0 to 1.0", "0.0");
    parser.add("--scale-full-probability", ArgType::Real, "Probability of scale full scheme operation, from 0.0 to 1.0", "0.1");
    parser.add("--scale-values", ArgType::String, "Scale coefficients", "1 -1 2 1/2");
}
