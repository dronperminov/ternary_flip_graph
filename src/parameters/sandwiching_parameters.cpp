#include "sandwiching_parameters.h"

void SandwichingParameters::parse(const ArgParser &parser) {
    probability = std::stod(parser["--sandwiching-probability"]);
    maxDenominator = std::stoi(parser["--sandwiching-max-denominator"]);
    minRows = std::stoi(parser["--sandwiching-min-rows"]);
    maxRows = std::stoi(parser["--sandwiching-max-rows"]);
    minColumns = std::stoi(parser["--sandwiching-min-columns"]);
    maxColumns = std::stoi(parser["--sandwiching-max-columns"]);
    values = parseFractions(parser["--sandwiching-values"]);

    if (minRows > maxRows)
        std::swap(minRows, maxRows);

    if (minColumns > maxColumns)
        std::swap(minColumns, maxColumns);

    nonZero.clear();

    for (const auto &fraction : values)
        if (fraction)
            nonZero.push_back(fraction);
}

std::ostream& operator<<(std::ostream& os, const SandwichingParameters &parameters) {
    if (parameters.probability == 0) {
        os << "Sandwiching operation: not used" << std::endl;
        return os;
    }

    os << "Sandwiching parameters:" << std::endl;
    os << "- probability: " << parameters.probability << std::endl;
    os << "- max denominator: " << parameters.maxDenominator << std::endl;
    os << "- rows: " << parameters.minRows << " .. " << parameters.maxRows << std::endl;
    os << "- columns: " << parameters.minColumns << " .. " << parameters.maxColumns << std::endl;

    os << "- values: [";
    for (size_t i = 0; i < parameters.values.size(); i++)
        os << (i > 0 ? ", " : "") << parameters.values[i].pretty();
    os << "]" << std::endl;
    return os;
}

void SandwichingParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0.5");
    parser.add("--sandwiching-max-denominator", ArgType::Natural, "Maximum denominator of sandwich matrices", "1");
    parser.add("--sandwiching-min-rows", ArgType::Natural, "Minimum number of changed rows", "1");
    parser.add("--sandwiching-max-rows", ArgType::Natural, "Maximum number of changed rows", "3");
    parser.add("--sandwiching-min-columns", ArgType::Natural, "Minimum number of changed columns", "2");
    parser.add("--sandwiching-max-columns", ArgType::Natural, "Maximum number of changed columns", "3");
    parser.add("--sandwiching-values", ArgType::String, "Sandwiching matrix values", "0 0 0 0 0 0 0 0 0 0 0 1 -1");
}
