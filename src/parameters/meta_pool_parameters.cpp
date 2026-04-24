#include "meta_pool_parameters.h"

void MetaPoolParameters::parse(const ArgParser &parser) {
    use = parser.isSet("--use-pool");
    size = parseNatural(parser["--pool-size"]);
    uniqueOnly = parser.isSet("--pool-unique-only");
    resume = parser.isSet("--resume");
    liftOnly = parser.isSet("--lift-only");
    alternativesProbability = std::stod(parser["--save-alternatives-probability"]);

    mergeMaxDiff = std::stoi(parser["--merge-max-diff"]);
    extendMaxDiff = std::stoi(parser["--extend-max-diff"]);
    projectMaxDiff = std::stoi(parser["--project-max-diff"]);
    productMaxDiff = std::stoi(parser["--product-max-diff"]);

    mergeProbability = std::stod(parser["--merge-probability"]);
    extendProbability = std::stod(parser["--extend-probability"]);
    projectProbability = std::stod(parser["--project-probability"]);
    productProbability = std::stod(parser["--product-probability"]);

    projectMinN1 = std::stoi(parser["--project-min-n1"]);
    projectMinN2 = std::stoi(parser["--project-min-n2"]);
    projectMinN3 = std::stoi(parser["--project-min-n3"]);

    selectRankScale = std::stod(parser["--select-rank-scale"]);
    metaRankScale = std::stod(parser["--meta-rank-scale"]);

    prioritiesPath = parser["--priorities-path"];
}

void MetaPoolParameters::writeJSON(std::ostream &os) const {
    os << "{";
    os << "\"size\": " << size << ", ";
    os << "\"unique_only\": " << (uniqueOnly ? "true" : "false") << ", ";
    os << "\"resume\": " << (resume ? "true" : "false") << ", ";
    os << "\"lift_only\": " << (liftOnly ? "true" : "false") << ", ";
    os << "\"save_alternatives_probability\": " << alternativesProbability << ", ";
    os << "\"merge_max_diff\": " << mergeMaxDiff << ", ";
    os << "\"extend_max_diff\": " << extendMaxDiff << ", ";
    os << "\"project_max_diff\": " << projectMaxDiff << ", ";
    os << "\"product_max_diff\": " << productMaxDiff << ", ";
    os << "\"merge_probability\": " << mergeProbability << ", ";
    os << "\"extend_probability\": " << extendProbability << ", ";
    os << "\"project_probability\": " << projectProbability << ", ";
    os << "\"product_probability\": " << productProbability << ", ";
    os << "\"project_min_n\": [" << projectMinN1 << ", " << projectMinN2 << ", " << projectMinN3 << "], ";
    os << "\"select_rank_scale\": " << selectRankScale << ", ";
    os << "\"meta_rank_scale\": " << metaRankScale << ", ";
    os << "\"priorities_path\": " << prioritiesPath;
    os << "}";
}

std::ostream& operator<<(std::ostream& os, const MetaPoolParameters &parameters) {
    if (parameters.use) {
        os << "Pool parameters:" << std::endl;
        os << "- size: " << parameters.size << std::endl;
        os << "- unique only: " << (parameters.uniqueOnly ? "yes" : "no") << std::endl;
        os << "- liftable only: " << (parameters.liftOnly ? "yes" : "no") << std::endl;
        os << "- resume: " << (parameters.resume ? "yes" : "no") << std::endl;
        os << "- save alternatives probability: " << parameters.alternativesProbability << std::endl;
        os << "- merge parameters (max diff: " << parameters.mergeMaxDiff << ", probability: " << parameters.mergeProbability << ")" << std::endl;
        os << "- extend parameters (max diff: " << parameters.extendMaxDiff << ", probability: " << parameters.extendProbability << ")" << std::endl;
        os << "- project parameters (max diff: " << parameters.projectMaxDiff << ", probability: " << parameters.projectProbability << ")" << std::endl;
        os << "- product parameters (max diff: " << parameters.productMaxDiff << ", probability: " << parameters.productProbability << ")" << std::endl;
        os << "- project min N: " << parameters.projectMinN1 << ", " << parameters.projectMinN2 << ", " << parameters.projectMinN3 << std::endl;
        os << "- select rank scale: " << parameters.selectRankScale << std::endl;
        os << "- meta rank scale: " << parameters.metaRankScale << std::endl;
        os << "- priorities path: " << parameters.prioritiesPath << std::endl;
    }

    return os;
}

void MetaPoolParameters::addToParser(ArgParser &parser, const std::string &sectionName) {
    parser.addSection(sectionName);
    parser.add("--use-pool", ArgType::Flag, "Use pool strategy");
    parser.add("--pool-size", ArgType::Natural, "Optimal size of pool", "1K");
    parser.add("--pool-unique-only", ArgType::Flag, "Save only unique schemes");
    parser.add("--resume", ArgType::Flag, "Read schemes from output directories as initial");
    parser.add("--lift-only", ArgType::Flag, "Save only schemes that can lift");
    parser.add("--save-alternatives-probability", ArgType::Real, "Save alternative schemes after runner end probability", "0.01");

    parser.add("--merge-max-diff", ArgType::UInt, "Max rank difference for merge operator", "4");
    parser.add("--extend-max-diff", ArgType::UInt, "Max rank difference for extend operator", "4");
    parser.add("--project-max-diff", ArgType::UInt, "Max rank difference for project operator", "4");
    parser.add("--product-max-diff", ArgType::UInt, "Max rank difference for product operator", "4");

    parser.add("--merge-probability", ArgType::Real, "Probability of merge operator", "0.25");
    parser.add("--extend-probability", ArgType::Real, "Probability of extend operator", "1.0");
    parser.add("--project-probability", ArgType::Real, "Probability of project operator", "0.8");
    parser.add("--product-probability", ArgType::Real, "Probability of product operator", "0.1");

    parser.add("--project-min-n1", ArgType::Natural, "Min first dimension for project", "2");
    parser.add("--project-min-n2", ArgType::Natural, "Min second dimension for project", "3");
    parser.add("--project-min-n3", ArgType::Natural, "Min third dimension for project", "4");

    parser.add("--select-rank-scale", ArgType::Real, "Scale for pool rank selection", "0.7");
    parser.add("--meta-rank-scale", ArgType::Real, "Scale for make meta operator", "0.5");

    parser.add("--priorities-path", ArgType::String, "Path to file with priorities", "priorities.txt");
}
