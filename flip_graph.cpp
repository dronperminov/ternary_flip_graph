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
#include "src/flip_graph.hpp"

int getMaxMatrixElements(const ArgParser &parser) {
    int n1, n2, n3;

    if (parser.isSet("--input-path")) {
        std::string path = parser["--input-path"];
        std::ifstream f(path);

        if (!f) {
            std::cerr << "Unable to open file \"" << path << "\"" << std::endl;
            return 0;
        }

        int count = 1;

        if (parser["--multiple"] == "true")
            f >> count;

        f >> n1 >> n2 >> n3;
        f.close();
    }
    else {
        n1 = std::stoi(parser["-n1"]);
        n2 = std::stoi(parser["-n2"]);
        n3 = std::stoi(parser["-n3"]);
    }

    int nn = std::max(n1 * n2, std::max(n2 * n3, n3 * n1));

    if (nn < 1) {
        std::cerr << "input matrix sizes is not valid: " << n1 << "x" << n2 << "x" << n3 << std::endl;
        return 0;
    }

    if (nn > 64) {
        std::cerr << "input matrix sizes too big: " << n1 << "x" << n2 << "x" << n3 << std::endl;
        return 0;
    }

    return nn;
}

template <template<typename> typename Scheme, typename T>
int runFlipGraph(const ArgParser &parser) {
    std::string outputPath = parser["--output-path"];

    std::string ring = parser["--ring"];
    int targetRank = std::stoi(parser["--target-rank"]);
    size_t flipIterations = parseNatural(parser["--flip-iterations"]);
    size_t minPlusIterations = parseNatural(parser["--min-plus-iterations"]);
    size_t maxPlusIterations = parseNatural(parser["--max-plus-iterations"]);
    size_t resetIterations = parseNatural(parser["--reset-iterations"]);
    int plusDiff = std::stoi(parser["--plus-diff"]);
    double sandwichingProbability = std::stod(parser["--sandwiching-probability"]);
    double reduceProbability = std::stod(parser["--reduce-probability"]);
    double copyBestProbability = std::stod(parser["--copy-best-probability"]);
    int maxImprovements = std::stoi(parser["--max-improvements"]);

    int count = std::stoi(parser["--count"]);
    int threads = std::stoi(parser["--threads"]);
    int topCount = std::stoi(parser["--top-count"]);
    int seed = std::stoi(parser["--seed"]);
    std::string format = parser["--format"];

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the flip graph algorithm:" << std::endl;
    if (parser.isSet("--input-path"))
        std::cout << "- input path: " << parser["--input-path"] << std::endl;
    else
        std::cout << "- dimension: " << parser["-n1"] << "x" << parser["-n2"] << "x" << parser["-n3"] << std::endl;

    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- target rank: " << targetRank << std::endl;
    std::cout << "- flip iterations: " << flipIterations << std::endl;
    std::cout << "- plus iterations: " << minPlusIterations << " .. " << maxPlusIterations << std::endl;
    std::cout << "- reset iterations: " << resetIterations << std::endl;
    std::cout << "- plus diff: " << plusDiff << std::endl;
    std::cout << "- sandwiching probability: " << sandwichingProbability << std::endl;
    std::cout << "- reduce probability: " << reduceProbability << std::endl;
    std::cout << "- copy best probability: " << copyBestProbability << std::endl;
    std::cout << "- max improvements: " << maxImprovements << std::endl;
    std::cout << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << std::endl;

    if (!makeDirectory(outputPath))
        return -1;

    FlipGraph<Scheme<T>> flipGraph(count, outputPath, threads, flipIterations, minPlusIterations, maxPlusIterations, resetIterations, plusDiff, sandwichingProbability, reduceProbability, copyBestProbability, seed, topCount, maxImprovements, format);

    bool valid;
    if (parser.isSet("--input-path")) {
        valid = flipGraph.initializeFromFile(parser["--input-path"], parser["--multiple"] == "true");
    }
    else {
        valid = flipGraph.initializeNaive(std::stoi(parser["-n1"]), std::stoi(parser["-n2"]), std::stoi(parser["-n3"]));
    }

    if (!valid)
        return -1;

    flipGraph.run(targetRank);
    return 0;
}

template <template<typename> typename Scheme>
int runFlipGraphSizes(const ArgParser &parser, int nn) {
    if (nn <= 16)
        return runFlipGraph<Scheme, uint16_t>(parser);

    if (nn <= 32)
        return runFlipGraph<Scheme, uint32_t>(parser);

    return runFlipGraph<Scheme, uint64_t>(parser);
}

int main(int argc, char **argv) {
    ArgParser parser("flip_graph", "Find fast matrix multiplication schemes using flip graph");

    parser.addSection("Matrix dimensions");
    parser.add("-n1", ArgType::Natural, "Number of rows in first matrix (A)");
    parser.add("-n2", ArgType::Natural, "Number of columns in A / rows in second matrix (B)");
    parser.add("-n3", ArgType::Natural, "Number of columns in second matrix (B)");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme(s)");
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for discovered schemes", "schemes");
    parser.add("--multiple", "-m", ArgType::Flag, "Read multiple schemes from file, with total count on first line");

    parser.addSection("Flip graph parameters");
    parser.addChoices("--ring", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2} or ZT - {-1, 0, 1}", {"ZT", "Z2", "Z3"}, "ZT");
    parser.add("--target-rank", ArgType::Natural, "Stop search when this rank is found, 0 searches for minimum", "0");
    parser.add("--flip-iterations", ArgType::Natural, "Flip iterations before reporting", "100K");
    parser.add("--min-plus-iterations", ArgType::Natural, "Minimum period for plus operator calls", "5K");
    parser.add("--max-plus-iterations", ArgType::Natural, "Maximum period for plus operator calls", "100K");
    parser.add("--reset-iterations", ArgType::Natural, "Total iterations before reset", "100M");
    parser.add("--plus-diff", ArgType::Natural, "Maximum rank difference for plus operations", "4");
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0");
    parser.add("--reduce-probability", ArgType::Real, "Probability of reduce operation, from 0.0 to 1.0", "0");
    parser.add("--copy-best-probability", ArgType::Real, "Probability to replace scheme with best scheme after improvement, from 0.0 to 1.0", "0.5");
    parser.add("--max-improvements", ArgType::Natural, "Maximum saved recent improvements for reset sampling", "10");

    parser.addSection("Run parameters");
    parser.add("--count", "-c", ArgType::Natural, "Number of parallel runners", "8");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.add("--top-count", ArgType::Natural, "Number of top schemes to report", "10");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.addChoices("--format", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    if (!parser.parse(argc, argv))
        return 0;

    if (!parser.isSet("--input-path") && (!parser.isSet("-n1") || !parser.isSet("-n2") || !parser.isSet("-n3"))) {
        std::cerr << "Must provide either dimension args (-n1 -n2 -n3) or an input file (-i)" << std::endl;
        return -1;
    }

    if (parser.isSet("--input-path") && (parser.isSet("-n1") || parser.isSet("-n2") || parser.isSet("-n3"))) {
        std::cerr << "Specify either dimension args (-n1 -n2 -n3) or an input file (-i), not both" << std::endl;
        return -1;
    }

    if (!parser.isSet("--input-path") && parser.isSet("--multiple")) {
        std::cerr << "--multiple flag requires an input file (-i), not dimension flags" << std::endl;
        return -1;
    }

    int nn = getMaxMatrixElements(parser);
    if (!nn)
        return -1;

    if (parser["--ring"] == "Z2")
        return runFlipGraphSizes<BinaryScheme>(parser, nn);

    if (parser["--ring"] == "Z3")
        return runFlipGraphSizes<Mod3Scheme>(parser, nn);

    return runFlipGraphSizes<TernaryScheme>(parser, nn);
}
