#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <omp.h>

#include "src/entities/arg_parser.h"
#include "src/entities/flip_graph.hpp"
#include "src/schemes/ternary_scheme.hpp"
#include "src/schemes/binary_scheme.hpp"

template <typename Scheme>
bool runFlipGraph(const ArgParser &parser, const std::string &ring) {
    int n1 = std::stoi(parser.get("-n1"));
    int n2 = std::stoi(parser.get("-n2"));
    int n3 = std::stoi(parser.get("-n3"));
    std::string inputPath = parser.get("-i");
    std::string outputPath = parser.get("-o");
    int count = std::stoi(parser.get("--count"));
    int targetRank = std::stoi(parser.get("--target-rank"));
    size_t flipIterations = std::stoul(parser.get("--flip-iterations"));
    size_t plusIterations = std::stoul(parser.get("--plus-iterations"));
    size_t resetIterations = std::stoul(parser.get("--reset-iterations"));
    int plusDiff = std::stoi(parser.get("--plus-diff"));
    double reduceProbability = std::stod(parser.get("--reduce-probability"));
    int seed = std::stoi(parser.get("--seed"));
    int threads = std::stoi(parser.get("--threads"));
    int topCount = std::stoi(parser.get("--top-count"));

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the ternary flip graph algorithm:" << std::endl;
    std::cout << "- dimension: " << n1 << "x" << n2 << "x" << n3 << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- target rank: " << targetRank << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus iterations: " << plusIterations << std::endl;
    std::cout << "- reset iterations: " << resetIterations << std::endl;
    std::cout << "- plus diff: " << plusDiff << std::endl;
    std::cout << "- reduce probability: " << reduceProbability << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    std::cout << std::endl;

    Scheme scheme;

    if (inputPath == "NULL") {
        if (!scheme.initializeNaive(n1, n2, n3))
            return false;
    }
    else if (!scheme.read(inputPath)) {
        return false;
    }

    if (!scheme.validate()) {
        std::cout << "error: initial scheme is invalid" << std::endl;
        return false;
    }

    std::cout << "Initial scheme correct, start!" << std::endl;
    std::cout << std::endl;

    FlipGraph<Scheme> flipGraph(count, outputPath, threads, flipIterations, plusIterations, resetIterations, plusDiff, reduceProbability, seed, topCount);
    flipGraph.run(scheme, targetRank);
    return true;
}

int main(int argc, char **argv) {
    ArgParser parser("ternary_flip_graph", "Find fast matrix multiplication using ternary flip graph");

    parser.add("-n1", ArgType::Natural, "INT", "number of first matrix rows");
    parser.add("-n2", ArgType::Natural, "INT", "number of first matrix columns and second matrix rows");
    parser.add("-n3", ArgType::Natural, "INT", "number of second matrix columns");
    parser.add("-i", ArgType::String, "PATH", "path to input file with scheme", "NULL");
    parser.add("-o", ArgType::String, "PATH", "path to save schemes", "schemes");
    parser.add("--count", ArgType::Natural, "INT", "number of parallel runners", "8");
    parser.add("--target-rank", ArgType::Natural, "INT", "rank for stop finding", "0");
    parser.add("--flip-iterations", ArgType::Natural, "INT", "number of iterations if flip call before reporting ", "10000");
    parser.add("--plus-iterations", ArgType::Natural, "INT", "number of iterations of plus call", "8000");
    parser.add("--reset-iterations", ArgType::Natural, "INT", "number of iterations of plus call", "100000000");
    parser.add("--plus-diff", ArgType::Natural, "INT", "maximum difference between scheme rank and best rank for plus call", "1");
    parser.add("--reduce-probability", ArgType::Real, "REAL", "probability of reduce call", "0");
    parser.add("--ring", ArgType::String, "Z2/ZT", "working ring (Z2 or ZT)", "ZT");
    parser.add("--seed", ArgType::Natural, "INT", "random seed", "0");
    parser.add("--threads", ArgType::Natural, "INT", "number of threads", std::to_string(omp_get_max_threads()));
    parser.add("--top-count", ArgType::Natural, "INT", "count of records for reporting", "10");

    if (!parser.parse(argc, argv))
        return 0;

    std::string ring = parser.get("--ring");

    if (ring == "Z2")
        return runFlipGraph<BinaryScheme<uint64_t>>(parser, ring);

    if (ring == "ZT")
        return runFlipGraph<TernaryScheme<uint64_t>>(parser, ring);

    std::cout << "error: invalid ring option: \"" << ring << "\"" << std::endl;
    return -1;
}
