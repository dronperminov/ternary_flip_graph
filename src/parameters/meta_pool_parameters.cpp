#include "meta_pool_parameters.h"

void MetaPoolParameters::parse(const ArgParser &parser) {
    use = parser.isSet("--use-pool");
    size = parseNatural(parser["--pool-size"]);
    uniqueOnly = parser.isSet("--pool-unique-only");
    alternativesProbability = std::stod(parser["--save-alternatives-probability"]);

    mergeMaxDiff = std::stoi(parser["--merge-max-diff"]);
    extendMaxDiff = std::stoi(parser["--extend-max-diff"]);
    projectMaxDiff = std::stoi(parser["--project-max-diff"]);
    productMaxDiff = std::stoi(parser["--product-max-diff"]);

    projectMinN1 = std::stoi(parser["--project-min-n1"]);
    projectMinN2 = std::stoi(parser["--project-min-n2"]);
    projectMinN3 = std::stoi(parser["--project-min-n3"]);

    selectRankScale = std::stod(parser["--select-rank-scale"]);
    metaRankScale = std::stod(parser["--meta-rank-scale"]);
}

void MetaPoolParameters::writeJSON(std::ostream &os) const {
    os << "{";
    os << "\"size\": " << size << ", ";
    os << "\"unique_only\": " << (uniqueOnly ? "true" : "false") << ", ";
    os << "\"save_alternatives_probability\": " << alternativesProbability << ", ";
    os << "\"merge_max_diff\": " << mergeMaxDiff << ", ";
    os << "\"extend_max_diff\": " << extendMaxDiff << ", ";
    os << "\"project_max_diff\": " << projectMaxDiff << ", ";
    os << "\"product_max_diff\": " << productMaxDiff << ", ";
    os << "\"project_min_n\": [" << projectMinN1 << ", " << projectMinN2 << ", " << projectMinN3 << "], ";
    os << "\"select_rank_scale\": " << selectRankScale << ", ";
    os << "\"meta_rank_scale\": " << metaRankScale << ", ";
    os << "}";
}

std::ostream& operator<<(std::ostream& os, const MetaPoolParameters &parameters) {
    if (parameters.use) {
        os << "Pool parameters:" << std::endl;
        os << "- size: " << parameters.size << std::endl;
        os << "- unique only: " << (parameters.uniqueOnly ? "yes" : "no") << std::endl;
        os << "- save alternatives probability: " << parameters.alternativesProbability << std::endl;
        os << "- merge max diff: " << parameters.mergeMaxDiff << std::endl;
        os << "- extend max diff: " << parameters.extendMaxDiff << std::endl;
        os << "- project max diff: " << parameters.projectMaxDiff << std::endl;
        os << "- product max diff: " << parameters.productMaxDiff << std::endl;
        os << "- project min N: " << parameters.projectMinN1 << ", " << parameters.projectMinN2 << ", " << parameters.projectMinN3 << std::endl;
        os << "- select rank scale: " << parameters.selectRankScale << std::endl;
        os << "- meta rank scale: " << parameters.metaRankScale << std::endl;
    }

    return os;
}

void MetaPoolParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--use-pool", ArgType::Flag, "Use pool strategy");
    parser.add("--pool-size", ArgType::Natural, "Optimal size of pool", "1K");
    parser.add("--pool-unique-only", ArgType::Flag, "Save only unique schemes");
    parser.add("--save-alternatives-probability", ArgType::Real, "Save alternative schemes after runner end probability", "0.1");

    parser.add("--merge-max-diff", ArgType::UInt, "Max rank difference for merge operator", "5");
    parser.add("--extend-max-diff", ArgType::UInt, "Max rank difference for extend operator", "10");
    parser.add("--project-max-diff", ArgType::UInt, "Max rank difference for project operator", "10");
    parser.add("--product-max-diff", ArgType::UInt, "Max rank difference for product operator", "5");

    parser.add("--project-min-n1", ArgType::Natural, "Min first dimension for project", "4");
    parser.add("--project-min-n2", ArgType::Natural, "Min second dimension for project", "5");
    parser.add("--project-min-n3", ArgType::Natural, "Min third dimension for project", "6");

    parser.add("--select-rank-scale", ArgType::Real, "Scale for pool rank selection", "0.7");
    parser.add("--meta-rank-scale", ArgType::Real, "Scale for make meta operator", "0.5");
}
