#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <random>
#include <ctime>
#include <unordered_set>
#include <omp.h>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/entities/sha1.h"
#include "src/schemes/binary_scheme.hpp"
#include "src/schemes/mod3_scheme.hpp"
#include "src/schemes/ternary_scheme.hpp"

template <typename Scheme>
std::string getSavePath(const Scheme &scheme, const std::string &outputPath, int version) {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2);
    SHA1 sha1;
    ss << "_m" << scheme.getRank();
    ss << "_" << sha1.get(scheme.getHash());
    ss << "_" << scheme.getRing();
    return ss.str();
}

template <typename Scheme>
std::string getHash(const Scheme &scheme, const std::string &unique) {
    if (unique == "structure")
        return scheme.getStructureHash();

    return scheme.getHash();
}

template <template<typename> typename Scheme, typename T>
int runFindAlternativeSchemes(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    std::string ring = parser["--ring"];
    std::string unique = parser["--unique-check"];
    double sandwichingProbability = std::stod(parser["--sandwiching-probability"]);
    double plusProbability = std::stod(parser["--plus-probability"]);
    int plusDiff = std::stoi(parser["--plus-diff"]);
    int threads = std::stoi(parser["--threads"]);

    size_t maxCount = parseNatural(parser["--max-count"]);
    int seed = std::stoi(parser["--seed"]);
    std::string format = parser["--format"];
    int maxMatrixElements = sizeof(T) * 8;

    if (seed == 0)
        seed = time(0);

    std::cout << "Start finding alternative schemes" << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- unique check: " << unique << std::endl;
    std::cout << std::endl;
    std::cout << "- sandwiching probability: " << sandwichingProbability << std::endl;
    std::cout << "- plus probability: " << plusProbability << std::endl;
    std::cout << "- plus diff: " << plusDiff << std::endl;
    std::cout << std::endl;
    std::cout << "- max count: " << maxCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << "- max matrix elements: " << maxMatrixElements << " (uint" << maxMatrixElements << "_t)" << std::endl;
    std::cout << std::endl << std::endl;

    Scheme<T> initialScheme;
    if (!initialScheme.read(inputPath, !parser.isSet("--no-verify")))
        return -1;

    if (!makeDirectory(outputPath))
        return -1;

    int schemeRank = initialScheme.getRank();
    int targetRank = parser.isSet("--target-rank") ? std::stoi(parser["--target-rank"]) : schemeRank;

    std::cout << "Readed scheme parameters:" << std::endl;
    std::cout << "- dimension: " << initialScheme.getDimension() << std::endl;
    std::cout << "- rank: " << initialScheme.getRank();

    if (targetRank != schemeRank)
        std::cout << " (target: " << targetRank << ")";

    std::cout << std::endl;
    std::cout << std::endl;

    std::vector<std::mt19937> generators = initRandomGenerators(seed, threads);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::unordered_set<std::string> hashes;
    hashes.insert(getHash(initialScheme, unique));

    std::vector<Scheme<T>> schemes(threads);
    for (int i = 0; i < threads; i++)
        schemes[i].copy(initialScheme);

    size_t count = 0;

    std::cout << "+-----------+-------------+------------+-----------------+" << std::endl;
    std::cout << "| iteration | alternative | complexity | available flips |" << std::endl;
    std::cout << "+-----------+-------------+------------+-----------------+" << std::endl;

    for (size_t iteration = 1; count < maxCount; iteration++) {
        std::vector<std::string> hashesThreads(threads, "");

        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < threads; i++) {
            Scheme<T>& scheme = schemes[i];
            std::mt19937& generator = generators[i];

            if (!scheme.tryFlip(generator) || (scheme.getRank() < targetRank + plusDiff && uniform(generator) < plusProbability))
                scheme.tryExpand(generator);

            if (uniform(generator) < sandwichingProbability)
                scheme.trySandwiching(generator);

            if (scheme.getRank() != targetRank)
                continue;

            hashesThreads[i] = getHash(scheme, unique);
        }

        for (int i = 0; i < threads; i++) {
            std::string hash = hashesThreads[i];
            if (hash.empty() || hashes.find(hash) != hashes.end())
                continue;

            count++;
            std::cout << "| ";
            std::cout << std::setw(9) << iteration << " | ";
            std::cout << std::setw(11) << count << " | ";
            std::cout << std::setw(10) << schemes[i].getComplexity() << " | ";
            std::cout << std::setw(15) << schemes[i].getAvailableFlips() << " |";
            std::cout << std::endl;

            std::string path = getSavePath(schemes[i], outputPath, count);
            if (format == "json")
                schemes[i].saveJson(path + ".json");
            else
                schemes[i].saveTxt(path + ".txt");

            hashes.insert(hash);
        }
    }

    std::cout << "+-----------+-------------+------------+-----------------+" << std::endl;
    return 0;
}

template <template<typename> typename Scheme>
int runFindAlternativeSchemesSizes(const ArgParser &parser) {
    int maxMatrixElements = getMaxMatrixElements(parser["--input-path"], false);
    if (maxMatrixElements < 0)
        return -1;

    if (maxMatrixElements <= 16)
        return runFindAlternativeSchemes<Scheme, uint16_t>(parser);

    if (maxMatrixElements <= 32)
        return runFindAlternativeSchemes<Scheme, uint32_t>(parser);

    if (maxMatrixElements <= 64)
        return runFindAlternativeSchemes<Scheme, uint64_t>(parser);

    if (maxMatrixElements <= 128)
        return runFindAlternativeSchemes<Scheme, __uint128_t>(parser);

    return runFindAlternativeSchemes<Scheme, uint256_t>(parser);
}

int main(int argc, char **argv) {
    ArgParser parser("find_alternative_schemes", "Find alternative fast matrix multiplication schemes using flip graph");
    parser.addChoices("--ring", "-r", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2} or ZT - {-1, 0, 1}", {"ZT", "Z2", "Z3"}, "ZT");
    parser.addChoices("--format", "-f", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "txt");
    parser.add("--max-count", "-n", ArgType::Natural, "Number of alternative schemes", "10K");
    parser.addChoices("--unique-check", "-u", ArgType::String, "Unique schemes comparator", {"coefficients", "structure"}, "coefficients");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for alternative schemes", "schemes");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Flip graph parameters");
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0.0");
    parser.add("--plus-probability", ArgType::Real, "Probability of plus operation, from 0.0 to 1.0", "0.001");
    parser.add("--plus-diff", ArgType::Natural, "Maximum rank difference for plus operations", "2");

    parser.addSection("Other parameters");
    parser.add("--target-rank", ArgType::Natural, "Rank of alternative schemes");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");

    if (!parser.parse(argc, argv))
        return 0;

    if (parser["--ring"] == "Z2")
        return runFindAlternativeSchemesSizes<BinaryScheme>(parser);

    if (parser["--ring"] == "Z3")
        return runFindAlternativeSchemesSizes<Mod3Scheme>(parser);

    return runFindAlternativeSchemesSizes<TernaryScheme>(parser);
}
