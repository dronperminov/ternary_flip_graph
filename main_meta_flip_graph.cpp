#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <omp.h>

#include "src/entities/arg_parser.h"
#include "src/schemes/ternary_scheme.hpp"
#include "src/schemes/binary_scheme.hpp"
#include "src/meta_flip_graph.hpp"

size_t parseNatural(std::string value) {
    size_t multiplier = 1;

    if (value.back() == 'K') {
        multiplier = 1000;
    }
    else if (value.back() == 'M') {
        multiplier = 1000000;
    }

    if (multiplier > 1)
        value.pop_back();

    return std::stoul(value) * multiplier;
}

template <template<typename> typename Scheme, typename T>
bool runMetaFlipGraph(const ArgParser &parser, int n1, int n2, int n3, const std::string inputPath, const std::string &ring) {
    std::string outputPath = parser.get("-o");
    int count = std::stoi(parser.get("--count"));
    size_t flipIterations = parseNatural(parser.get("--flip-iterations"));
    size_t minPlusIterations = parseNatural(parser.get("--min-plus-iterations"));
    size_t maxPlusIterations = parseNatural(parser.get("--max-plus-iterations"));
    size_t resetIterations = parseNatural(parser.get("--reset-iterations"));
    int plusDiff = std::stoi(parser.get("--plus-diff"));
    double reduceProbability = std::stod(parser.get("--reduce-probability"));
    double resizeProbability = std::stod(parser.get("--resize-probability"));
    int seed = std::stoi(parser.get("--seed"));
    int threads = std::stoi(parser.get("--threads"));
    int topCount = std::stoi(parser.get("--top-count"));

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the ternary meta flip graph algorithm:" << std::endl;
    if (inputPath == "NULL")
        std::cout << "- dimension: " << n1 << "x" << n2 << "x" << n3 << std::endl;
    else
        std::cout << "- input path: " << inputPath << std::endl;

    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus iterations: " << minPlusIterations << " .. " << maxPlusIterations << std::endl;
    std::cout << "- reset iterations: " << resetIterations << std::endl;
    std::cout << "- plus diff: " << plusDiff << std::endl;
    std::cout << "- reduce probability: " << reduceProbability << std::endl;
    std::cout << "- resize probability: " << resizeProbability << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << std::endl;

    MetaFlipGraph<Scheme<T>> metaFlipGraph(count, outputPath, threads, flipIterations, minPlusIterations, maxPlusIterations, resetIterations, plusDiff, reduceProbability, resizeProbability, seed, topCount);

    if (ring == "ZT" || ring == "ternary")
        metaFlipGraph.initializeBestTernaryRanks();
    else if (ring == "Z2" || ring == "binary")
        metaFlipGraph.initializeBestBinaryRanks();

    if (inputPath == "NULL") {
        metaFlipGraph.initializeNaive(n1, n2, n3);
    }
    else if (!metaFlipGraph.initializeFromFile(inputPath))
        return false;

    metaFlipGraph.run();
    return true;
}

int main(int argc, char **argv) {
    ArgParser parser("ternary_meta_flip_graph", "Find fast matrix multiplication schemes using ternary flip graph");

    parser.add("-n1", ArgType::Natural, "INT", "number of rows in first matrix (A)", "0");
    parser.add("-n2", ArgType::Natural, "INT", "number of columns in A / rows in second matrix (B)", "0");
    parser.add("-n3", ArgType::Natural, "INT", "number of columns in second matrix (B)", "0");

    parser.add("-i", ArgType::String, "PATH", "path to input file with initial scheme", "NULL");
    parser.add("-o", ArgType::String, "PATH", "output directory for discovered schemes", "schemes");

    parser.add("--flip-iterations", ArgType::Natural, "INT", "flip iterations before reporting ", "100K");
    parser.add("--min-plus-iterations", ArgType::Natural, "INT", "minimum plus iterations per phase", "5K");
    parser.add("--max-plus-iterations", ArgType::Natural, "INT", "maximum plus iterations per phase", "100K");
    parser.add("--reset-iterations", ArgType::Natural, "INT", "total iterations before reset", "100M");
    parser.add("--plus-diff", ArgType::Natural, "INT", "maximum rank difference for plus operations", "4");
    parser.add("--reduce-probability", ArgType::Real, "REAL", "probability of reduce operation (0.0 to 1.0)", "0");
    parser.add("--resize-probability", ArgType::Real, "REAL", "probability of resize operation (0.0 to 1.0)", "0");
    parser.add("--ring", ArgType::String, "Z2/ZT", "coefficient ring: Z2 ({0, 1}) or ZT ({-1, 0, 1})", "ZT");

    parser.add("--count", ArgType::Natural, "INT", "number of parallel runners", "8");
    parser.add("--threads", ArgType::Natural, "INT", "number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.add("--top-count", ArgType::Natural, "INT", "number of top schemes to report", "10");
    parser.add("--seed", ArgType::Natural, "INT", "random seed (0 = time-based)", "0");

    if (!parser.parse(argc, argv))
        return 0;

    int n1 = std::stoi(parser.get("-n1"));
    int n2 = std::stoi(parser.get("-n2"));
    int n3 = std::stoi(parser.get("-n3"));
    std::string inputPath = parser.get("-i");

    if (inputPath == "NULL" && (n1 == 0 || n2 == 0 || n3 == 0)) {
        std::cout << "error: must provide either dimension flags (-n1 -n2 -n3) or an input file (-i)" << std::endl;
        return -1;
    }

    if (inputPath != "NULL" && (n1 != 0 || n2 != 0 || n3 != 0)) {
        std::cout << "error: specify either dimension flags (-n1 -n2 -n3) or an input file (-i), not both" << std::endl;
        return -1;
    }

    std::string ring = parser.get("--ring");

    if (ring == "Z2" || ring == "binary")
        return runMetaFlipGraph<BinaryScheme, uint64_t>(parser, n1, n2, n3, inputPath, ring);

    if (ring == "ZT" || ring == "ternary")
        return runMetaFlipGraph<TernaryScheme, uint64_t>(parser, n1, n2, n3, inputPath, ring);

    std::cout << "error: invalid ring option: \"" << ring << "\"" << std::endl;
    return -1;
}
