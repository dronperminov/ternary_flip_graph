#include "plus_parameters.h"

void PlusParameters::parse(const ArgParser &parser) {
    probability = std::stod(parser["--plus-probability"]);
    iterations = std::stoi(parser["--plus-iterations"]);
    values = parseFractions(parser["--split-values"]);
}

std::ostream& operator<<(std::ostream& os, const PlusParameters &parameters) {
    if (parameters.probability == 0) {
        os << "Plus operation: not used" << std::endl;
        return os;
    }

    os << "Plus parameters:" << std::endl;
    os << "- probability: " << parameters.probability << std::endl;
    os << "- iterations: " << parameters.iterations << std::endl;
    os << "- split values: [";
    for (size_t i = 0; i < parameters.values.size(); i++)
        os << (i > 0 ? ", " : "") << parameters.values[i].pretty();
    os << "]" << std::endl;
    return os;
}

void PlusParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--plus-probability", ArgType::Real, "Probability of plus operation, from 0.0 to 1.0", "0.01");
    parser.add("--plus-iterations", ArgType::Real, "Max iterations of flips after plus operation", "1000");
    parser.add("--split-values", ArgType::String, "Split values", "0 1 -1");
}
