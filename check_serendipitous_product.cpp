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
#include "src/entities/sha1.h"
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

std::unordered_map<std::string, int> getKnownRanks(const std::string &ring) {
    std::unordered_map<std::string, int> dimension2ranks = KNOWN_RANKS.at(ring);
    std::unordered_map<std::string, int> knownRanks;

    for (int i = 1; i <= 16; i++) {
        for (int j = 1; j <= 16; j++) {
            for (int k = 1; k <= 16; k++) {
                std::string key = getDimension(i, j, k, true);
                std::string dimension = getDimension(i, j, k);

                if (dimension2ranks.find(key) == dimension2ranks.end())
                    knownRanks[dimension] = i * j * k;
                else
                    knownRanks[dimension] = dimension2ranks.at(key);
            }
        }
    }

    return knownRanks;
}

bool checkSerendipitousProduct(const FractionalScheme &scheme, std::mt19937 &generator, std::unordered_map<std::string, int> &knownRanks, int iterations) {
    FlipStructureOptimizer optimizer = scheme.getStructureOptimizer();
    int n[3] = {
        scheme.getDimension(0),
        scheme.getDimension(1),
        scheme.getDimension(2)
    };

    bool found = false;

    for (int n1 = 1; n1 <= 16 && n[0] * n1 <= 16; n1++) {
        for (int n2 = 1; n2 <= 16 && n[1] * n2 <= 16; n2++) {
            for (int n3 = 1; n3 <= 16 && n[2] * n3 <= 16; n3++) {
                if (n1 == 1 && n2 == 1 && n3 == 1)
                    continue;

                int d[3] = {n1, n2, n3};
                int serendipitousRank = optimizer.getSerendipitousRank(generator, d, knownRanks, iterations);
                if (serendipitousRank == -1)
                    continue;

                std::string dimension = getDimension(n1 * n[0], n2 * n[1], n3 * n[2], true);
                int knownRank = knownRanks.at(dimension);

                if (serendipitousRank < knownRank) {
                    std::cout << "Serendipitous for " << dimension << ": " << serendipitousRank << " < " << knownRank << std::endl;
                    knownRanks[dimension] = serendipitousRank;
                    found = true;
                }
            }
        }
    }

    return found;
}

int runCheckSerendipitousProduct(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];
    bool shuffle = parser.isSet("--shuffle-paths");
    std::string ring = parser["--ring"];
    std::string format = parser["--format"];

    int threads = std::stoi(parser["--threads"]);
    bool verify = !parser.isSet("--no-verify");
    int iterations = std::stoi(parser["--iterations"]);

    int seed = std::stoi(parser["--seed"]);
    if (seed == 0)
        seed = time(0);

    if (!std::filesystem::is_directory(inputPath) && !(std::filesystem::is_regular_file(inputPath) && endsWith(inputPath, ".txt"))) {
        std::cout << "Input path must be a directory or .txt file with paths" << std::endl;
        return -1;
    }

    if (!makeDirectory(outputPath))
        return -1;

    std::vector<std::mt19937> generators = initRandomGenerators(seed, threads);

    std::cout << "Parsed parameters of the check_serendipitous_product tool:" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- shuffle paths: " << (shuffle ? "yes" : "no") << std::endl;
    std::cout << "- ring: " << ring << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- iterations: " << iterations << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << std::endl;

    std::vector<std::string> paths = getSchemePaths(inputPath, shuffle, generators[0]);
    if (paths.empty()) {
        std::cout << "There are no scheme files" << std::endl;
        return 0;
    }

    std::cout << "Found " << paths.size() << " files" << std::endl;
    int digits = digitsCount(paths.size());

    std::unordered_map<std::string, int> knownRanks = getKnownRanks(ring);

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < paths.size(); i++) {
        std::string path = paths[i];

        std::stringstream ss;
        ss << std::setw(digits) << (i + 1) << " / " << paths.size() << ". " << path << std::endl;
        std::cout << ss.str();

        FractionalScheme scheme;
        if (!scheme.read(path, verify, !endsWith(path, "Q.txt")))
            continue;

        if (checkSerendipitousProduct(scheme, generators[omp_get_thread_num()], knownRanks, iterations))
            scheme.save(outputPath + "/" + scheme.getFilename(format));
    }

    return 0;
}

int main(int argc, char *argv[]) {
    ArgParser parser("check_serendipitous_product", "Check existing schemes for serendipitous product");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.addChoices("--ring", "-r", ArgType::String, "Ring for improvements checking", {"ZT", "Z", "Q"}, "Q");
    parser.addChoices("--format", "-f", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "txt");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input directory with schemes or .txt file with paths for check", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for found schemes", "schemes");
    parser.add("--shuffle-paths", ArgType::Flag, "Shuffle readed paths");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Other parameters");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.add("--iterations", ArgType::Natural, "Iterations of random checking", "25");

    if (!parser.parse(argc, argv))
        return 0;

    return runCheckSerendipitousProduct(parser);
}
