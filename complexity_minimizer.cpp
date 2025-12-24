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
#include "src/schemes/binary_scheme.hpp"
#include "src/complexity_minimizer.hpp"

template <template<typename> typename Scheme, typename T>
int runComplexityMinimizer(const ArgParser &parser) {
    std::string inputPath = parser.get("-i");
    std::string outputPath = parser.get("-o");

    size_t flipIterations = parseNatural(parser.get("--flip-iterations"));
    double plusProbability = std::stod(parser.get("--plus-probability"));
    std::string ring = parser.get("--ring");

    int count = std::stoi(parser.get("--count"));
    int threads = std::stoi(parser.get("--threads"));
    int topCount = std::stoi(parser.get("--top-count"));
    int seed = std::stoi(parser.get("--seed"));
    int maxNoImprovements = std::stoi(parser.get("--max-no-improvements"));
    std::string format = parser.get("--format");

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the complexity minimizer algorithm:" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus probability: " << plusProbability << std::endl;
    std::cout << "- ring: " << ring << std::endl;
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

    parser.add("-i", ArgType::String, "PATH", "path to input file with initial scheme");
    parser.add("-o", ArgType::String, "PATH", "output directory for discovered schemes", "schemes");

    parser.add("--flip-iterations", ArgType::Natural, "INT", "flip iterations before reporting ", "100K");
    parser.add("--plus-probability", ArgType::Real, "REAL", "probability of plus operation (0.0 to 1.0)", "0");
    parser.add("--ring", ArgType::String, "Z2/ZT", "coefficient ring: Z2 ({0, 1}) or ZT ({-1, 0, 1})", "ZT");

    parser.add("--count", ArgType::Natural, "INT", "number of parallel runners", "8");
    parser.add("--threads", ArgType::Natural, "INT", "number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.add("--top-count", ArgType::Natural, "INT", "number of top schemes to report", "10");
    parser.add("--seed", ArgType::Natural, "INT", "random seed (0 = time-based)", "0");
    parser.add("--max-no-improvements", ArgType::Natural, "INT", "maximum iterations without complexity improvement before termination", "3");
    parser.add("--format", ArgType::String, "txt/json", "output format for saved schemes: txt or json", "json");

    if (!parser.parse(argc, argv))
        return 0;

    std::string ring = parser.get("--ring");

    if (ring == "Z2" || ring == "binary")
        return runComplexityMinimizer<BinaryScheme, uint64_t>(parser);

    if (ring == "ZT" || ring == "ternary")
        return runComplexityMinimizer<TernaryScheme, uint64_t>(parser);

    std::cout << "error: invalid ring option: \"" << ring << "\"" << std::endl;
    return -1;
}
