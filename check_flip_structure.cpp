#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <omp.h>
#include <filesystem>
#include <random>

#include "src/utils.h"
#include "src/known_ranks.h"
#include "src/entities/arg_parser.h"
#include "src/entities/flip_structure_optimizer.h"
#include "src/entities/buffer_writer.h"
#include "src/schemes/fractional_scheme.h"

std::vector<std::string> getSchemePaths(const std::string &inputPath, bool shuffle, std::mt19937 &generator) {
    std::vector<std::string> paths;
    std::vector<std::string> extensions = {"ZT.txt", "Z.txt", "Q.txt"};

    if (std::filesystem::is_directory(inputPath)) {
        std::cout << "Start reading files from directory \"" << inputPath << "\"" << std::endl;
        paths = getSchemePathsFromDirectoryRecursive(inputPath, extensions);
    }
    else {
        std::cout << "Start reading paths file \"" << inputPath << "\"" << std::endl;
        paths = getSchemePathsFromFile(inputPath, extensions);
    }

    if (shuffle)
        std::shuffle(paths.begin(), paths.end(), generator);

    return paths;
}

std::string getLogPath(const std::string &outputPath, int thread) {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << "structure_log" << std::setw(2) << std::setfill('0') << thread << ".jsonl";
    return ss.str();
}

int runCheckFlipStructure(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];
    bool shuffle = parser.isSet("--shuffle-paths");
    bool onlyOptimal = parser.isSet("--only-optimal");

    int threads = std::stoi(parser["--threads"]);
    bool verify = !parser.isSet("--no-verify");
    int iterations = std::stoi(parser["--iterations"]);
    double eps = std::stod(parser["--eps"]);

    int seed = std::stoi(parser["--seed"]);
    if (seed == 0)
        seed = time(0);

    if (!std::filesystem::is_directory(inputPath) && !(std::filesystem::is_regular_file(inputPath) && endsWith(inputPath, ".txt"))) {
        std::cout << "Input path must be a directory or .txt file with paths" << std::endl;
        return -1;
    }

    if (!makeDirectory(outputPath))
        return -1;

    std::cout << "Parsed parameters of the check_structure tool:" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- shuffle paths: " << (shuffle ? "yes" : "no") << std::endl;
    std::cout << "- only optimal: " << (onlyOptimal ? "yes" : "no") << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- iterations: " << iterations << std::endl;
    std::cout << "- eps: " << eps << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << std::endl;

    std::vector<std::mt19937> generators = initRandomGenerators(seed, threads);
    std::vector<std::string> paths = getSchemePaths(inputPath, shuffle, generators[0]);
    if (paths.empty()) {
        std::cout << "There are no scheme files" << std::endl;
        return 0;
    }

    std::cout << "Found " << paths.size() << " files" << std::endl;

    std::unordered_map<std::string, int> dimension2rank = KNOWN_RANKS.at("Q");
    std::vector<BufferWriter> writers;

    for (int i = 0; i < threads; i++)
        writers.emplace_back(BufferWriter(getLogPath(outputPath, i), 512));

    std::cout << "+-------------+-----------+------+-------------------+-------------------+" << std::endl;
    std::cout << "| path number | dimension | rank |       omega       |  structure omega  |" << std::endl;
    std::cout << "+-------------+-----------+------+-------------------+-------------------+" << std::endl;

    #pragma omp parallel for schedule(dynamic, std::max(1, std::min(64, (int)paths.size() / threads))) num_threads(threads)
    for (size_t i = 0; i < paths.size(); i++) {
        std::string path = paths[i];
        int thread = omp_get_thread_num();

        FractionalScheme scheme;
        if (!scheme.read(path, verify, !endsWith(path, "Q.txt")))
            continue;

        std::string dimension = scheme.getDimension();
        if (dimension2rank.find(dimension) == dimension2rank.end() || (onlyOptimal && scheme.getRank() > dimension2rank.at(dimension)))
            continue;

        FlipStructureOptimizer optimizer = scheme.getFullStructureOptimizer();
        FlipStructure structure = optimizer.optimize(generators[thread], iterations, eps);

        std::stringstream ss;
        ss << "| ";
        ss << std::setw(11) << (i + 1) << " | ";
        ss << std::setw(9) << dimension << " | ";
        ss << std::setw(4) << scheme.getRank() << " | ";
        ss << std::setw(17) << std::setprecision(15) << scheme.getOmega() << " | ";
        ss << std::setw(17) << std::setprecision(15) << structure.omega << " |";
        ss << std::endl;
        std::cout << ss.str();

        std::stringstream line;
        line << "{";
        line << "\"path\": \"" << path << "\", ";
        line << "\"dimension\": [" << scheme.getDimension(0) << ", " << scheme.getDimension(1) << ", " << scheme.getDimension(2) << "], ";
        line << "\"rank\": " << scheme.getRank() << ", ";
        line << "\"omega\": " << std::setprecision(15) << scheme.getOmega() << ", ";
        line << "\"structure_omega\": " << std::setprecision(15) << structure.omega << ", ";
        line << "\"structure\": " << structure;
        line << "}";

        writers[thread].add(line.str());        
    }

    std::cout << "+-------------+-----------+------+-------------------+-------------------+" << std::endl;

    return 0;
}

int main(int argc, char *argv[]) {
    ArgParser parser("check_serendipitous_product", "Check existing schemes for serendipitous product");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.add("--only-optimal", ArgType::Flag, "Check only schemes with optimal rank");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input directory with schemes or .txt file with paths for check", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for found schemes", "schemes/structure");
    parser.add("--shuffle-paths", ArgType::Flag, "Shuffle readed paths");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Other parameters");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.add("--iterations", ArgType::Natural, "Iterations of random checking", "100");
    parser.add("--eps", ArgType::Real, "Precision of omega", "1e-10");

    if (!parser.parse(argc, argv))
        return 0;

    return runCheckFlipStructure(parser);
}
