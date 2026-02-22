#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <omp.h>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/schemes/ternary_scheme.hpp"
#include "src/schemes/mod3_scheme.hpp"
#include "src/schemes/binary_scheme.hpp"
#include "src/scheme_optimizer.hpp"

template <template<typename> typename Scheme, typename T>
int runSchemeOptimizer(const ArgParser &parser) {
    std::string ring = parser["--ring"];
    int count = std::stoi(parser["--count"]);
    int threads = std::stoi(parser["--threads"]);
    std::string format = parser["--format"];

    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    size_t flipIterations = parseNatural(parser["--flip-iterations"]);
    double plusProbability = std::stod(parser["--plus-probability"]);
    int plusDiff = std::stoi(parser["--plus-diff"]);

    int topCount = std::stoi(parser["--top-count"]);
    int seed = std::stoi(parser["--seed"]);
    std::string metric = parser["--metric"];
    bool maximize = parser.isSet("--maximize");
    double copyBestProbability = std::stod(parser["--copy-best-probability"]);
    int maxNoImprovements = std::stoi(parser["--max-no-improvements"]);
    int maxMatrixElements = sizeof(T) * 8;

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the " << metric << " " << (maximize ? "maximizer" : "minimizer") << " algorithm:" << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus probability: " << plusProbability << std::endl;
    std::cout << "- plus diff: " << plusDiff << std::endl;
    std::cout << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- copy best probability: " << copyBestProbability << std::endl;
    std::cout << "- max no improvements: " << maxNoImprovements << std::endl;
    std::cout << "- max matrix elements: " << maxMatrixElements << " (uint" << maxMatrixElements << "_t)" << std::endl;
    std::cout << std::endl;

    SchemeOptimizer<Scheme<T>> optimizer(count, outputPath, threads, flipIterations, plusProbability, plusDiff, seed, copyBestProbability, metric, maximize, topCount, format);

    if (!optimizer.initializeFromFile(inputPath, parser.isSet("--multiple"), !parser.isSet("--no-verify")))
        return -1;

    if (!makeDirectory(outputPath))
        return -1;

    optimizer.run(maxNoImprovements);
    return true;
}

template <template<typename> typename Scheme>
int runSchemeOptimizerSizes(const ArgParser &parser) {
    int maxMatrixElements = getMaxMatrixElements(parser["--input-path"], parser.isSet("--multiple"));
    if (maxMatrixElements < 0)
        return -1;

    if (maxMatrixElements <= 16)
        return runSchemeOptimizer<Scheme, uint16_t>(parser);

    if (maxMatrixElements <= 32)
        return runSchemeOptimizer<Scheme, uint32_t>(parser);

    if (maxMatrixElements <= 64)
        return runSchemeOptimizer<Scheme, uint64_t>(parser);

    return runSchemeOptimizer<Scheme, __uint128_t>(parser);
}

int main(int argc, char **argv) {
    ArgParser parser("scheme_optimizer", "Optimize fast matrix multiplication schemes for naive complexity or potential flips count using flip graph");

    parser.addChoices("--ring", "-r", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2} or ZT - {-1, 0, 1}", {"ZT", "Z2", "Z3"}, "ZT");
    parser.add("--count", "-c", ArgType::Natural, "Number of parallel runners", "8");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.addChoices("--format", "-f", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "txt");
    parser.addChoices("--metric", ArgType::String, "Metric for optimization", {"complexity", "flips"}, "complexity");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme(s)", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for optimized schemes", "schemes");
    parser.add("--multiple", "-m", ArgType::Flag, "Read multiple schemes from file, with total count on first line");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Random walk parameters");
    parser.add("--flip-iterations", ArgType::Natural, "Flip iterations before reporting", "100K");
    parser.add("--plus-probability", ArgType::Real, "Probability of plus operation, from 0.0 to 1.0", "0.01");
    parser.add("--plus-diff", ArgType::Natural, "Maximum rank difference for plus operations", "2");

    parser.addSection("Other parameters");
    parser.add("--top-count", ArgType::Natural, "Number of top schemes to report", "10");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.add("--maximize", ArgType::Flag, "Maximize instead of minimizing");
    parser.add("--copy-best-probability", ArgType::Real, "Probability to replace scheme with best scheme after improvement, from 0.0 to 1.0", "0.5");
    parser.add("--max-no-improvements", ArgType::Natural, "Maximum iterations without metric improvement before termination", "3");

    if (!parser.parse(argc, argv))
        return 0;

    if (parser["--ring"] == "Z2")
        return runSchemeOptimizerSizes<BinaryScheme>(parser);

    if (parser["--ring"] == "Z3")
        return runSchemeOptimizerSizes<Mod3Scheme>(parser);

    return runSchemeOptimizerSizes<TernaryScheme>(parser);
}
