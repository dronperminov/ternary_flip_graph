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

std::string getWriterPath(const std::string &outputPath, int thread, int digits) {
    std::stringstream ss;
    ss << outputPath << "/" << "analyzed_thread" << std::setw(digits) << std::setfill('0') << thread << ".jsonl";
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const std::vector<Flip> &flips) {
    os << "[";

    for (int p = 0; p < 3; p++) {
        os << (p > 0 ? ", " : "") << "[";

        bool printed = false;
        for (const Flip& flip : flips) {
            if (flip.p != p)
                continue;

            os << (printed ? ", " : "") << "[" << flip.i << ", " << flip.j << "]";
            printed = true;
        }

        os << "]";
    }

    os << "]";
    return os;
}

int runAnalyzeSchemes(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    bool shuffle = parser.isSet("--shuffle-paths");
    int threads = std::stoi(parser["--threads"]);
    bool verify = !parser.isSet("--no-verify");
    int seed = time(0);

    bool saveComplexity = parser.isSet("--save-complexity");
    bool saveRing = parser.isSet("--save-ring");
    bool saveBuds = parser.isSet("--save-buds");
    bool saveCoefficientSet = parser.isSet("--save-coefficient-set");
    bool saveTypeInvariant = parser.isSet("--save-type-invariant");
    bool saveBudsInvariant = parser.isSet("--save-buds-invariant");
    bool saveSignCanonized = parser.isSet("--save-sign-canonized");
    bool saveAll = parser.isSet("--save-all");

    if (!std::filesystem::is_directory(inputPath) && !(std::filesystem::is_regular_file(inputPath) && endsWith(inputPath, ".txt"))) {
        std::cout << "Input path must be a directory or .txt file with paths" << std::endl;
        return -1;
    }

    if (!makeDirectory(outputPath))
        return -1;

    std::cout << "Parsed parameters of the analyze_schemes tool:" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- shuffle paths: " << (shuffle ? "yes" : "no") << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << std::endl;
    std::cout << "Analyzed parameters:" << std::endl;
    std::cout << "- complexity: " << (saveComplexity || saveAll ? "yes" : "no") << std::endl;
    std::cout << "- ring: " << (saveRing || saveAll ? "yes" : "no") << std::endl;
    std::cout << "- buds: " << (saveBuds || saveAll ? "yes" : "no") << std::endl;
    std::cout << "- coefficient set: " << (saveCoefficientSet || saveAll ? "yes" : "no") << std::endl;
    std::cout << "- type invariant: " << (saveTypeInvariant || saveAll ? "yes" : "no") << std::endl;
    std::cout << "- buds invariant: " << (saveBudsInvariant || saveAll ? "yes" : "no") << std::endl;

    std::mt19937 generator(seed);
    std::vector<std::string> paths = getSchemePaths(inputPath, shuffle, generator);
    if (paths.empty()) {
        std::cout << "There are no scheme files" << std::endl;
        return 0;
    }

    std::cout << "Found " << paths.size() << " files" << std::endl;

    int digits = digitsCount(threads);
    std::vector<BufferWriter> writers;
    for (int i = 0; i < threads; i++)
        writers.emplace_back(BufferWriter(getWriterPath(outputPath, i, digits), 512));

    std::cout << "+-------------+-----------+------+-------------------+" << std::endl;
    std::cout << "| path number | dimension | rank |       omega       |" << std::endl;
    std::cout << "+-------------+-----------+------+-------------------+" << std::endl;

    #pragma omp parallel for schedule(dynamic, std::max(1, std::min(64, (int)paths.size() / threads))) num_threads(threads)
    for (size_t i = 0; i < paths.size(); i++) {
        std::string path = paths[i];
        int thread = omp_get_thread_num();

        FractionalScheme scheme;
        if (!scheme.read(path, verify, !endsWith(path, "Q.txt")))
            continue;

        std::stringstream ss;
        ss << "| ";
        ss << std::setw(11) << (i + 1) << " | ";
        ss << std::setw(9) << scheme.getDimension() << " | ";
        ss << std::setw(4) << scheme.getRank() << " | ";
        ss << std::setw(17) << std::setprecision(15) << scheme.getOmega() << " | ";
        ss << std::endl;
        std::cout << ss.str();

        std::replace(path.begin(), path.end(), '\\', '/');

        std::stringstream line;
        line << "{";
        line << "\"path\": \"" << path << "\"";
        line << ", \"dimension\": [" << scheme.getDimension(0) << ", " << scheme.getDimension(1) << ", " << scheme.getDimension(2) << "]";
        line << ", \"rank\": " << scheme.getRank() << "";
        line << ", \"omega\": " << std::setprecision(15) << scheme.getOmega();

        if (saveComplexity || saveAll)
            line << ", \"complexity\": " << scheme.getComplexity();

        if (saveRing || saveAll)
            line << ", \"ring\": \"" << scheme.getRing() << "\"";

        if (saveBuds || saveBudsInvariant || saveAll) {
            FlipStructureOptimizer optimizer = scheme.getFullStructureOptimizer();

            if (saveBuds || saveAll)
                line << ", \"buds\": " << optimizer.getFlips();

            if (saveBudsInvariant || saveAll)
                line << ", \"buds_invariant\": \"" << optimizer.getBudsInvariant() << "\"";
        }

        if (saveCoefficientSet || saveAll)
            line << ", \"coefficient_set\": \"" << scheme.getUniqueValues() << "\"";

        if (saveTypeInvariant || saveAll)
            line << ", \"type_invariant\": \"" << scheme.getTypeInvariant() << "\"";

        if (saveSignCanonized || saveAll)
            line << ", \"sign_canonized\": " << (scheme.isSignCanonized() ? "true" : "false");

        line << "}";

        writers[thread].add(line.str());        
    }

    std::cout << "+-------------+-----------+------+-------------------+" << std::endl;

    return 0;
}

int main(int argc, char *argv[]) {
    ArgParser parser("analyze_schemes", "Get information about schemes");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input directory with schemes or .txt file with paths for check", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for results", "schemes/analyzed");
    parser.add("--shuffle-paths", ArgType::Flag, "Shuffle readed paths");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Analyzed parameters");
    parser.add("--save-complexity", "-c", ArgType::Flag, "Save naive addition complexity");
    parser.add("--save-ring", "-r", ArgType::Flag, "Save naive addition complexity");
    parser.add("--save-buds", "-b", ArgType::Flag, "Save buds indices for U, V and W parts");
    parser.add("--save-coefficient-set", "-cs", ArgType::Flag, "Save unique coefficient set");
    parser.add("--save-type-invariant", "-ti", ArgType::Flag, "Save type invariant");
    parser.add("--save-buds-invariant", "-bi", ArgType::Flag, "Save buds invariant");
    parser.add("--save-sign-canonized", "-sc", ArgType::Flag, "Save sign canonized check");
    parser.add("--save-all", "-a", ArgType::Flag, "Save all available parameters");

    if (!parser.parse(argc, argv))
        return 0;

    return runAnalyzeSchemes(parser);
}
