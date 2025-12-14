#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>

#include "src/schemes/ternary_scheme.hpp"
#include "src/entities/arg_parser.h"


std::string getSavePath(const TernaryScheme &scheme, int iteration, const std::string path) {
    std::stringstream ss;
    ss << path << "/";
    ss << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2);
    ss << "_m" << scheme.getRank();
    ss << "_c" << scheme.getComplexity();
    ss << "_iteration" << iteration;
    ss << ".json";

    return ss.str();
}


std::string prettyInt(int value) {
    std::stringstream ss;

    if (value < 1000)
        ss << value;
    else if (value < 1000000)
        ss << std::setprecision(2) << std::fixed << (value / 1000.0) << "K";
    else
        ss << std::setprecision(2) << std::fixed << (value / 1000000.0) << "M";

    return ss.str();
}

int main(int argc, char **argv) {
    ArgParser parser("ternary_flip_graph", "Find fast matrix multiplication using ternary flip graph");

    parser.add("-n1", ArgType::Natural, "INT", "number of first matrix rows");
    parser.add("-n2", ArgType::Natural, "INT", "number of first matrix columns and second matrix rows");
    parser.add("-n3", ArgType::Natural, "INT", "number of second matrix columns");
    parser.add("-i", ArgType::String, "PATH", "path to input file with scheme", "NULL");
    parser.add("-o", ArgType::String, "PATH", "path to save schemes", "schemes");
    parser.add("--target-rank", ArgType::Natural, "INT", "rank for stop finding", "0");
    parser.add("--plus-iterations", ArgType::Natural, "INT", "number of iterations of plus call", "8000");
    parser.add("--reduce-probability", ArgType::Real, "REAL", "probability of reduce call", "0");
    parser.add("--log-period", ArgType::Natural, "INT", "period of report printing", "10000");
    parser.add("--seed", ArgType::Natural, "INT", "random seed", "0");

    if (!parser.parse(argc, argv))
        return 0;

    int n1 = std::stoi(parser.get("-n1"));
    int n2 = std::stoi(parser.get("-n2"));
    int n3 = std::stoi(parser.get("-n3"));
    std::string inputPath = parser.get("-i");
    std::string outputPath = parser.get("-o");
    int targetRank = std::stoi(parser.get("--target-rank"));
    int plusIterations = std::stoi(parser.get("--plus-iterations"));
    double reduceProbability = std::stod(parser.get("--reduce-probability"));
    int logPeriod = std::stoi(parser.get("--log-period"));
    int seed = std::stoi(parser.get("--seed"));

    if (seed == 0)
        seed = time(0);

    TernaryScheme scheme;

    if (inputPath == "NULL") {
        scheme.initializeNaive(n1, n2, n3);
    }
    else if (!scheme.read(inputPath)) {
        return -1;
    }

    std::cout << "Start ternary flip graph algorithm" << std::endl;

    std::cout << "- dimension: " << n1 << "x" << n2 << "x" << n3 << std::endl;
    std::cout << "- plus iterations: " << plusIterations << std::endl;
    std::cout << "- reduce probability: " << reduceProbability << std::endl;
    std::cout << "- target rank: " << targetRank << std::endl;
    std::cout << "- log period: " << logPeriod << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << std::endl;

    if (!scheme.validate()) {
        std::cout << "error: initial scheme is invalid" << std::endl;
        return -1;
    }

    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    int bestRank = scheme.getRank();
    int flips = 0;

    for (size_t iteration = 0; 1; iteration++) {
        int prevRank = scheme.getRank();

        if (!scheme.tryFlip(generator)) {
            if (scheme.tryPlus(generator))
                flips = 0;

            continue;
        }

        if (uniform(generator) < reduceProbability && scheme.tryReduce())
            flips = 0;

        int rank = scheme.getRank();
        if (rank < prevRank)
            flips = 0;

        flips++;

        if (rank < bestRank) {
            if (!scheme.validate()) {
                std::cout << prettyInt(iteration) << ". Invalid scheme!" << std::endl;
                break;
            }

            std::string path = getSavePath(scheme, iteration, outputPath);
            scheme.save(path);
            std::cout << prettyInt(iteration) << ". Rank was improved from " << bestRank << " to " << rank << std::endl;
            bestRank = rank;
        }

        if (rank <= targetRank)
            break;

        if (iteration % logPeriod == 0)
            std::cout << prettyInt(iteration) << ". Flips: " << prettyInt(flips) << ", rank: " << rank << " (best: " << bestRank << ")" << std::endl;

        if (flips >= plusIterations && rank < bestRank + 1 && scheme.tryPlus(generator))
            flips = 0;
    }

    return 0;
}