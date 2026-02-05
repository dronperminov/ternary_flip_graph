#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <omp.h>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/entities/flip_parameters.h"
#include "src/entities/meta_parameters.h"
#include "src/schemes/ternary_scheme.hpp"
#include "src/schemes/mod3_scheme.hpp"
#include "src/schemes/binary_scheme.hpp"
#include "src/meta_flip_graph.hpp"

template <template<typename> typename Scheme, typename T>
int runMetaFlipGraph(const ArgParser &parser) {
    std::string ring = parser["--ring"];
    int count = std::stoi(parser["--count"]);
    int threads = std::stoi(parser["--threads"]);
    std::string format = parser["--format"];

    std::string outputPath = parser["--output-path"];

    FlipParameters flipParameters;
    flipParameters.parse(parser);

    MetaParameters metaParameters;
    metaParameters.parse(parser);

    int seed = std::stoi(parser["--seed"]);
    int topCount = std::stoi(parser["--top-count"]);
    std::string improveRing = parser["--improve-ring"];
    int maxMatrixElements = sizeof(T) * 8;

    if (seed == 0)
        seed = time(0);

    std::cout << "Parsed parameters of the meta flip graph algorithm:" << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- format: " << format << std::endl;

    if (parser.isSet("--input-path"))
        std::cout << "- input path: " << parser["--input-path"] << std::endl;
    else
        std::cout << "- dimension: " << parser["-n1"] << "x" << parser["-n2"] << "x" << parser["-n3"] << std::endl;

    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << flipParameters << std::endl;
    std::cout << metaParameters << std::endl;

    std::cout << "Other parameters:" << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- top count: " << topCount << std::endl;
    if (improveRing != "")
        std::cout << "- improve ring: " << improveRing << std::endl;
    std::cout << "- max matrix elements: " << maxMatrixElements << " (uint" << maxMatrixElements << "_t)" << std::endl;
    std::cout << std::endl;

    if (!makeDirectory(outputPath))
        return -1;

    MetaFlipGraph<Scheme<T>> metaFlipGraph(count, outputPath, threads, flipParameters, metaParameters, seed, topCount, format);
    metaFlipGraph.initializeKnownRanks(improveRing);

    bool valid;
    if (parser.isSet("--input-path")) {
        valid = metaFlipGraph.initializeFromFile(parser["--input-path"], parser.isSet("--multiple"));
    }
    else {
        valid = metaFlipGraph.initializeNaive(std::stoi(parser["-n1"]), std::stoi(parser["-n2"]), std::stoi(parser["-n3"]));
    }

    if (!valid)
        return -1;

    metaFlipGraph.run();
    return 0;
}

template <template<typename> typename Scheme>
int runMetaFlipGraphSizes(const ArgParser &parser) {
    int maxMatrixElements = std::stoi(parser["--int-width"]);

    if (parser.isSet("-n1") || parser.isSet("-n2") || parser.isSet("-n3")) {
        int n1 = std::stoi(parser["-n1"]);
        int n2 = std::stoi(parser["-n2"]);
        int n3 = std::stoi(parser["-n3"]);

        maxMatrixElements = std::max(n1 * n2, std::max(n2 * n3, n3 * n1));
    }

    if (maxMatrixElements <= 16)
        return runMetaFlipGraph<Scheme, uint16_t>(parser);

    if (maxMatrixElements <= 32)
        return runMetaFlipGraph<Scheme, uint32_t>(parser);

    if (maxMatrixElements <= 64)
        return runMetaFlipGraph<Scheme, uint64_t>(parser);

    return runMetaFlipGraph<Scheme, __uint128_t>(parser);
}

int main(int argc, char **argv) {
    ArgParser parser("meta_flip_graph", "Find fast matrix multiplication schemes using meta flip graph");

    parser.addChoices("--ring", "-r", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2} or ZT - {-1, 0, 1}", {"ZT", "Z2", "Z3"}, "ZT");
    parser.add("--count", "-c", ArgType::Natural, "Number of parallel runners", "8");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.addChoices("--format", "-f", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "txt");

    parser.addSection("Matrix dimensions (only for naive initialization)");
    parser.add("-n1", ArgType::Natural, "Number of rows in first matrix (A)");
    parser.add("-n2", ArgType::Natural, "Number of columns in A / rows in second matrix (B)");
    parser.add("-n3", ArgType::Natural, "Number of columns in second matrix (B)");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme(s)");
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for discovered schemes", "schemes");
    parser.add("--multiple", "-m", ArgType::Flag, "Read multiple schemes from file, with total count on first line");

    parser.addSection("Random walk parameters");
    parser.add("--flip-iterations", ArgType::Natural, "Flip iterations before reporting", "1M");
    parser.add("--min-plus-iterations", ArgType::Natural, "Minimum period for plus operator calls", "5K");
    parser.add("--max-plus-iterations", ArgType::Natural, "Maximum period for plus operator calls", "100K");
    parser.add("--reset-iterations", ArgType::Natural, "Total iterations before reset", "10B");
    parser.add("--plus-diff", ArgType::Natural, "Maximum rank difference for plus operations", "4");
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0");
    parser.add("--reduce-probability", ArgType::Real, "Probability of reduce operation, from 0.0 to 1.0", "0");

    parser.addSection("Meta operations parameters");
    parser.add("--meta-probability", ArgType::Real, "Probability of call meta operations, from 0.0 to 1.0", "0");
    parser.addChoices("--meta-strategy", ArgType::String, "Strategy of meta operations", {"default", "proj", "ext"}, "default");
    parser.add("--meta-min-dimension", ArgType::Natural, "Min dimension for project meta operation", "2");
    parser.add("--meta-max-dimension", ArgType::Natural, "Max dimension for merge/extend meta operations", "16");
    parser.add("--meta-max-rank", ArgType::Natural, "Max rank for merge/extend meta operations", "350");
    parser.add("--meta-max-rank-diff", ArgType::Natural, "Max rank difference for reset to initial", "10");

    parser.addSection("Other parameters");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.add("--top-count", ArgType::Natural, "Number of top schemes to report", "10");
    parser.addChoices("--improve-ring", ArgType::String, "Only save schemes that improve known rank for this ring (saves all by default)", {"Z2", "ZT", "Q", ""}, "");
    parser.addChoices("--int-width", ArgType::String, "Integer bit width (16/32/64/128), determines maximum matrix elements", {"16", "32", "64", "128"}, "64");

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

    if (parser["--ring"] == "Z2")
        return runMetaFlipGraphSizes<BinaryScheme>(parser);

    if (parser["--ring"] == "Z3")
        return runMetaFlipGraphSizes<Mod3Scheme>(parser);

    return runMetaFlipGraphSizes<TernaryScheme>(parser);
}
