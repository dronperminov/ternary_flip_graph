#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <random>
#include <ctime>
#include <unordered_set>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/schemes/binary_scheme.hpp"
#include "src/schemes/mod3_scheme.hpp"
#include "src/schemes/ternary_scheme.hpp"

template <template<typename> typename Scheme, typename T>
std::string getSavePath(const Scheme<T> &scheme, const std::string &outputPath, int version) {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2);
    ss << "_m" << scheme.getRank();
    ss << "_v" << std::setw(6) << std::setfill('0') << version;
    ss << "_" << scheme.getRing();
    return ss.str();
}

template <template<typename> typename Scheme, typename T>
int runFindAlternativeSchemes(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    std::string ring = parser["--ring"];
    double sandwichingProbability = std::stod(parser["--sandwiching-probability"]);
    double plusProbability = std::stod(parser["--plus-probability"]);
    int plusDiff = std::stoi(parser["--plus-diff"]);

    size_t maxCount = parseNatural(parser["--max-count"]);
    int seed = std::stoi(parser["--seed"]);
    std::string format = parser["--format"];

    if (seed == 0)
        seed = time(0);

    std::cout << "Start finding alternative schemes" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << "- sandwiching probability: " << sandwichingProbability << std::endl;
    std::cout << "- plus probability: " << plusProbability << std::endl;
    std::cout << "- plus diff: " << plusDiff << std::endl;
    std::cout << std::endl;
    std::cout << "- max count: " << maxCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << std::endl << std::endl;

    Scheme<T> scheme;
    if (!scheme.read(inputPath))
        return -1;

    int schemeRank = scheme.getRank();
    int targetRank = parser.isSet("--target-rank") ? std::stoi(parser["--target-rank"]) : schemeRank;

    std::cout << "Readed scheme parameters:" << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- dimension: " << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2) << std::endl;
    std::cout << "- rank: " << scheme.getRank();

    if (targetRank != schemeRank)
        std::cout << " (target: " << targetRank << ")";

    std::cout << std::endl;
    std::cout << std::endl;

    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::unordered_set<std::string> hashes;
    hashes.insert(scheme.getHash());

    size_t count = 0;

    std::cout << "+-----------+-------------+------------+-----------------+" << std::endl;
    std::cout << "| iteration | alternative | complexity | available flips |" << std::endl;
    std::cout << "+-----------+-------------+------------+-----------------+" << std::endl;

    for (size_t iteration = 1; count < maxCount; iteration++) {
        if (!scheme.tryFlip(generator) || (scheme.getRank() < targetRank + plusDiff && uniform(generator) < plusProbability))
            scheme.tryExpand(generator);

        if (uniform(generator) < sandwichingProbability)
            scheme.trySandwiching(generator);

        if (scheme.getRank() != targetRank)
            continue;

        std::string hash = scheme.getHash();
        if (hashes.find(hash) != hashes.end())
            continue;

        count++;
        std::cout << "| ";
        std::cout << std::setw(9) << iteration << " | ";
        std::cout << std::setw(11) << count << " | ";
        std::cout << std::setw(10) << scheme.getComplexity() << " | ";
        std::cout << std::setw(15) << scheme.getAvailableFlips() << " |";
        std::cout << std::endl;

        std::string path = getSavePath(scheme, outputPath, count);
        if (format == "json")
            scheme.saveJson(path + ".json");
        else
            scheme.saveTxt(path + ".txt");

        hashes.insert(hash);
    }

    std::cout << "+-----------+-------------+------------+" << std::endl;
    return 0;
}

int main(int argc, char **argv) {
    ArgParser parser("find_alternative_schemes", "Find alternative fast matrix multiplication schemes using flip graph");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for alternative schemes", "schemes");

    parser.addSection("Find algorithm parameters");
    parser.addChoices("--ring", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2} or ZT - {-1, 0, 1}", {"ZT", "Z2", "Z3"}, "ZT");
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0.5");
    parser.add("--plus-probability", ArgType::Real, "Probability of plus operation, from 0.0 to 1.0", "0.2");
    parser.add("--plus-diff", ArgType::Natural, "Maximum rank difference for plus operations", "2");
    parser.add("--target-rank", ArgType::Natural, "Rank of alternative schemes");

    parser.addSection("Run parameters");
    parser.add("--max-count", ArgType::Natural, "Number of alternative schemes", "10000");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.addChoices("--format", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    if (!parser.parse(argc, argv))
        return 0;

    if (parser["--ring"] == "Z2")
        return runFindAlternativeSchemes<BinaryScheme, uint64_t>(parser);

    if (parser["--ring"] == "Z3")
        return runFindAlternativeSchemes<Mod3Scheme, uint64_t>(parser);

    return runFindAlternativeSchemes<TernaryScheme, uint64_t>(parser);
}
