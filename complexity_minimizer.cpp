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
#include "src/complexity_minimizer.hpp"

template <template<typename> typename Scheme, typename T>
int runComplexityMinimizer(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    std::string ring = parser["--ring"];
    size_t flipIterations = parseNatural(parser["--flip-iterations"]);
    double plusProbability = std::stod(parser["--plus-probability"]);

    int count = std::stoi(parser["--count"]);
    int threads = std::stoi(parser["--threads"]);
    int topCount = std::stoi(parser["--top-count"]);
    int seed = std::stoi(parser["--seed"]);
    int maxNoImprovements = std::stoi(parser["--max-no-improvements"]);
    std::string format = parser["--format"];

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the complexity minimizer algorithm:" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus probability: " << plusProbability << std::endl;
    std::cout << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- max no improvements: " << maxNoImprovements << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << std::endl;

    ComplexityMinimizer<Scheme<T>> minimizer(count, outputPath, threads, flipIterations, plusProbability, seed, topCount, format);

    if (!minimizer.initializeFromFile(inputPath))
        return -1;

    minimizer.run(maxNoImprovements);
    return true;
}

int main(int argc, char **argv) {
    ArgParser parser("complexity_minimizer", "Find fast matrix multiplication scheme with lowest naive complexity using flip graph");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial schemes", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for minimized schemes", "schemes");

    parser.addSection("Complexity minimizer parameters");
    parser.addChoices("--ring", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2} or ZT - {-1, 0, 1}", {"ZT", "Z2", "Z3"}, "ZT");
    parser.add("--flip-iterations", ArgType::Natural, "Flip iterations before reporting ", "100K");
    parser.add("--plus-probability", ArgType::Real, "Probability of plus operation, from 0.0 to 1.0", "0");

    parser.addSection("Run parameters");
    parser.add("--count", "-c", ArgType::Natural, "Number of parallel runners", "8");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.add("--top-count", ArgType::Natural, "Number of top schemes to report", "10");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.add("--max-no-improvements", ArgType::Natural, "Maximum iterations without complexity improvement before termination", "3");
    parser.addChoices("--format", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    if (!parser.parse(argc, argv))
        return 0;

    if (parser["--ring"] == "Z2")
        return runComplexityMinimizer<BinaryScheme, uint64_t>(parser);

    if (parser["--ring"] == "Z3")
        return runComplexityMinimizer<Mod3Scheme, uint64_t>(parser);

    return runComplexityMinimizer<TernaryScheme, uint64_t>(parser);
}
