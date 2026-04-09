#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <sstream>
#include <omp.h>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/parameters/sandwiching_parameters.h"
#include "src/parameters/sandwich_flip_parameters.h"
#include "src/parameters/scale_parameters.h"
#include "src/parameters/plus_parameters.h"
#include "src/sandwich_flip_optimizer.h"

bool runSandwichFlip(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];
    int count = std::stoi(parser["--count"]);
    int threads = std::min(count, std::stoi(parser["--threads"]));

    SandwichFlipParameters sandwichFlipParameters;
    sandwichFlipParameters.parse(parser);

    SandwichingParameters sandwichingParameters;
    sandwichingParameters.parse(parser);

    ScaleParameters scaleParameters;
    scaleParameters.parse(parser);

    PlusParameters plusParameters;
    plusParameters.parse(parser);

    size_t maxImprovements = parseNatural(parser["--max-improvements"]);
    size_t maxNoImprovements = parseNatural(parser["--max-no-improvements"]);
    int seed = std::stoi(parser["--seed"]);
    if (seed == 0)
        seed = time(0);

    std::string format = parser["--format"];

    std::cout << "Parsed parameters of the sandwich-flip tool:" << std::endl;
    std::cout << "- count: " << count << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;

    std::cout << sandwichFlipParameters << std::endl;
    std::cout << sandwichingParameters << std::endl;
    std::cout << scaleParameters << std::endl;
    std::cout << plusParameters << std::endl;

    std::cout << "Other parameters:" << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << "- max no improvements: " << prettyInt(maxNoImprovements) << std::endl;
    std::cout << "- max improvements: " << maxImprovements << std::endl;
    std::cout << std::endl;

    double sumP = sandwichingParameters.probability + scaleParameters.probability;
    if (sumP > 1) {
        std::cout << "Invalid probabilities: sandwiching + scale > 1" << std::endl;
        return false;
    }

    SandwichFlipOptimizer optimizer(count, outputPath, threads, sandwichFlipParameters, sandwichingParameters, scaleParameters, plusParameters, seed, maxImprovements, format);
    if (!optimizer.initializeFromFile(inputPath, !parser.isSet("--no-verify"), parser.isSet("--integer")))
        return false;

    optimizer.run(maxNoImprovements);
    return true;
}

int main(int argc, char *argv[]) {
    ArgParser parser("sandwich_flip", "Try to optimize scheme using random flips, scales and sandwiching");
    parser.add("--count", "-c", ArgType::Natural, "Number of parallel runners", "16");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for improved schemes", "schemes");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");
    parser.add("--integer", ArgType::Flag, "Read scheme as integer");

    SandwichFlipParameters::addToParser(parser, "Run parameters");
    SandwichingParameters::addToParser(parser, "Sandwiching parameters");
    ScaleParameters::addToParser(parser, "Scale parameters");
    PlusParameters::addToParser(parser, "Plus parameters");

    parser.addSection("Other parameters");
    parser.add("--max-improvements", ArgType::Natural, "Number of last improved schemes", "8");
    parser.add("--max-no-improvements", ArgType::Natural, "Maximum number of iterations without improvement", "10M");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.addChoices("--format", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    if (!parser.parse(argc, argv))
        return 0;

    if (!runSandwichFlip(parser))
        return -1;

    return 0;
}
