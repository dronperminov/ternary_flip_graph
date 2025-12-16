#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <omp.h>

#include "src/entities/arg_parser.h"
#include "src/entities/flip_graph.h"

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
    parser.add("--reduce-probability", ArgType::Real, "REAL", "probability of reduce call", "0");
    parser.add("--seed", ArgType::Natural, "INT", "random seed", "0");
    parser.add("--threads", ArgType::Natural, "INT", "number of threads", std::to_string(omp_get_max_threads()));
    parser.add("--top-count", ArgType::Natural, "INT", "count of records for reporting", "10");

    if (!parser.parse(argc, argv))
        return 0;

    int n1 = std::stoi(parser.get("-n1"));
    int n2 = std::stoi(parser.get("-n2"));
    int n3 = std::stoi(parser.get("-n3"));
    std::string inputPath = parser.get("-i");
    std::string outputPath = parser.get("-o");
    int count = std::stoi(parser.get("--count"));
    int targetRank = std::stoi(parser.get("--target-rank"));
    int flipIterations = std::stoi(parser.get("--flip-iterations"));
    int plusIterations = std::stoi(parser.get("--plus-iterations"));
    double reduceProbability = std::stod(parser.get("--reduce-probability"));
    int seed = std::stoi(parser.get("--seed"));
    int threads = std::stoi(parser.get("--threads"));
    int topCount = std::stoi(parser.get("--top-count"));

    if (seed == 0)
        seed = time(0);

    TernaryScheme scheme;

    if (inputPath == "NULL") {
        if (!scheme.initializeNaive(n1, n2, n3))
            return -1;
    }
    else if (!scheme.read(inputPath)) {
        return -1;
    }

    std::cout << "Start ternary flip graph algorithm" << std::endl;
    std::cout << "- dimension: " << n1 << "x" << n2 << "x" << n3 << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- target rank: " << targetRank << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus iterations: " << plusIterations << std::endl;
    std::cout << "- reduce probability: " << reduceProbability << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    std::cout << std::endl;

    if (!scheme.validate()) {
        std::cout << "error: initial scheme is invalid" << std::endl;
        return -1;
    }

    FlipGraph flipGraph(count, outputPath, threads, flipIterations, plusIterations, reduceProbability, seed, topCount);
    flipGraph.run(scheme, targetRank);
    return 0;
}
