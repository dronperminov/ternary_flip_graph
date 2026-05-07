#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include <omp.h>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/entities/sha1.h"
#include "src/schemes/binary_scheme.hpp"
#include "src/schemes/mod3_scheme.hpp"
#include "src/schemes/fractional_scheme.h"
#include "src/lift/binary_lifter.h"
#include "src/lift/mod3_lifter.h"

std::vector<std::string> getSchemePaths(const std::string &inputPath) {
    if (std::filesystem::is_regular_file(inputPath))
        return {inputPath};

    std::cout << "Start reading files from directory \"" << inputPath << "\"" << std::endl;
    std::vector<std::string> paths = getSchemePathsFromDirectory(inputPath, {".txt"});
    if (paths.empty())
        std::cout << "Directory is empty: nothing to lift" << std::endl;

    return paths;
}

template <template<typename> typename Scheme, typename T>
int runLiftSchemes(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    std::string ring = parser["--ring"];
    int steps = std::stoi(parser["--steps"]);
    bool canonize = parser.isSet("--canonize");
    bool fixFractions = parser.isSet("--fix-fractions");

    int threads = std::stoi(parser["--threads"]);
    std::string format = parser["--format"];
    int maxMatrixElements = sizeof(T) * 8;

    if (!makeDirectory(outputPath))
        return -1;

    std::cout << "Lift schemes from " << ring << " field to general" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- steps: " << steps << std::endl;
    std::cout << "- canonize: " << (canonize ? "yes" : "no") << std::endl;
    std::cout << "- fix fractions: " << (fixFractions ? "yes" : "no") << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << "- max matrix elements: " << maxMatrixElements << " (uint" << maxMatrixElements << "_t)" << std::endl;
    std::cout << std::endl << std::endl;

    std::vector<std::string> paths = getSchemePaths(inputPath);
    if (paths.empty())
        return 0;

    std::cout << "Start lift " << paths.size() << " schemes from \"" << inputPath << "\"" << std::endl;
    std::cout << std::endl;

    std::cout << "+--------+-----------+------+----------------------------+-------+--------------+" << std::endl;
    std::cout << "| scheme | dimension | rank |           status           | steps | elapsed time |" << std::endl;
    std::cout << "+--------+-----------+------+----------------------------+-------+--------------+" << std::endl;

    std::vector<double> elapsedTimes(paths.size(), 0);
    auto startTime = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < paths.size(); i++) {
        std::string status;
        int step = 0;

        auto t1 = std::chrono::high_resolution_clock::now();

        Scheme<T> scheme;
        if (!scheme.read(paths[i], !parser.isSet("--no-verify"))) {
            status = "invalid scheme";
            steps = -1;
        }
        else {
            FractionalScheme liftedScheme;

            bool reconstructed = scheme.reconstruct(liftedScheme) && liftedScheme.validateParallel();

            if (!reconstructed) {
                auto lifter = scheme.toLift();

                while (step < steps && !reconstructed && lifter.lift()) {
                    reconstructed = lifter.reconstruct(liftedScheme) && liftedScheme.validateParallel();
                    step++;
                }
            }

            if (reconstructed) {
                if (fixFractions)
                    liftedScheme.fixFractions();

                if (canonize)
                    liftedScheme.canonize();

                status = "reconstructed in " + liftedScheme.getRing();
                std::string path = outputPath + "/" + liftedScheme.getFilename(format);

                if (format == "txt")
                    liftedScheme.saveTxt(path);
                else
                    liftedScheme.saveJson(path);
            }
            else if (step == steps) {
                status = "no rational reconstruction";
            }
            else {
                status = "lifting failed";
            }
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes[i] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0;

        std::stringstream ss;
        ss << "| " << std::setw(6) << (i + 1) << " | ";
        ss << std::setw(9) << scheme.getDimension() << " | ";
        ss << std::setw(4) << scheme.getRank() << " | ";
        ss << std::setw(26) << status << " | ";
        ss << std::setw(5) << step << " | ";
        ss << std::setw(12) << prettyTime(elapsedTimes[i]) << " |" << std::endl;
        std::cout << ss.str();
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    double elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 1000.0;
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    std::cout << "+--------+-----------+------+----------------------------+-------+--------------+" << std::endl;
    std::cout << "- elapsed time (total / mean): " << prettyTime(elapsedTime) << " / " << prettyTime(meanTime) << std::endl;
    return 0;
}

template <template<typename> typename Scheme>
int runLiftSchemesSizes(const ArgParser &parser) {
    int maxMatrixElements = parser["--int-width"] == "auto" ? getMaxMatrixElements(parser["--input-path"], false) : std::stoi(parser["--int-width"]);
    if (maxMatrixElements < 0)
        return -1;

    if (maxMatrixElements <= 16)
        return runLiftSchemes<Scheme, uint16_t>(parser);

    if (maxMatrixElements <= 32)
        return runLiftSchemes<Scheme, uint32_t>(parser);

    if (maxMatrixElements <= 64)
        return runLiftSchemes<Scheme, uint64_t>(parser);

    if (maxMatrixElements <= 128)
        return runLiftSchemes<Scheme, __uint128_t>(parser);

    return runLiftSchemes<Scheme, uint256_t>(parser);
}

int main(int argc, char *argv[]) {
    ArgParser parser("lift", "Lift schemes from Z2/Z3 field to general");
    parser.addChoices("--ring", "-r", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2}", {"Z2", "Z3"}, "", true);
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.addChoices("--format", "-f", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with scheme or directory with schemes", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for lifted schemes", "schemes/lifted");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Lifting parameters");
    parser.add("--steps", "-k", ArgType::Natural, "Number of Hensel lifting steps", "10");
    parser.add("--canonize", "-c", ArgType::Flag, "Canonize reconstructed schemes");
    parser.add("--fix-fractions", ArgType::Flag, "Try to rescale fractions to integers");

    parser.addSection("Other parameters");
    parser.addChoices("--int-width", ArgType::String, "Integer bit width (16/32/64/128/256), determines maximum matrix elements", {"16", "32", "64", "128", "256", "auto"}, "auto");

    if (!parser.parse(argc, argv))
        return 0;

    if (parser["--ring"] == "Z2")
        return runLiftSchemesSizes<BinaryScheme>(parser);

    if (parser["--ring"] == "Z3")
        return runLiftSchemesSizes<Mod3Scheme>(parser);

    return 0;
}
