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
#include "src/schemes/binary_scheme.hpp"
#include "src/schemes/mod3_scheme.hpp"
#include "src/schemes/fractional_scheme.h"
#include "src/lift/binary_lifter.h"
#include "src/lift/mod3_lifter.h"

template <typename Scheme>
bool readSchemes(const std::string &inputPath, std::vector<Scheme> &schemes, bool multiple, bool checkCorrectness) {
    std::ifstream f(inputPath);
    if (!f) {
        std::cerr << "Unable to open file \"" << inputPath << "\"" << std::endl;
        return false;
    }

    int count = 1;
    if (multiple)
        f >> count;

    std::cout << "Start reading " << count << " schemes" << std::endl;
    schemes.resize(count);

    bool valid = true;
    for (int i = 0; i < count && valid; i++)
        valid = schemes[i].read(f, checkCorrectness);

    f.close();
    return valid;
}

std::string getSavePath(const FractionalScheme &scheme, int index, const std::string &outputPath, const std::string &format) {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_" << scheme.getRing();
    ss << "_v" << index;
    ss << "." << format;
    return ss.str();
}

template <typename Scheme>
int runLiftSchemes(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    std::string ring = parser["--ring"];
    int steps = std::stoi(parser["--steps"]);
    bool canonize = parser.isSet("--canonize");

    int threads = std::stoi(parser["--threads"]);
    std::string format = parser["--format"];

    if (!makeDirectory(outputPath))
        return -1;

    std::cout << "Lift schemes from " << ring << " field to general" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- steps: " << steps << std::endl;
    std::cout << "- canonize: " << (canonize ? "yes" : "no") << std::endl;
    std::cout << "- threads: " << threads << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << std::endl << std::endl;

    std::vector<Scheme> schemes;
    if (!readSchemes(inputPath, schemes, parser.isSet("--multiple"), !parser.isSet("--no-verify")))
        return -1;

    std::cout << "Successfully read " << schemes.size() << " schemes from \"" << inputPath << "\"" << std::endl;
    std::cout << std::endl;

    std::cout << "+--------+-----------+------+----------------------------+-------+--------------+" << std::endl;
    std::cout << "| scheme | dimension | rank |           status           | steps | elapsed time |" << std::endl;
    std::cout << "+--------+-----------+------+----------------------------+-------+--------------+" << std::endl;

    std::vector<double> elapsedTimes(schemes.size(), 0);
    auto startTime = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < schemes.size(); i++) {
        auto t1 = std::chrono::high_resolution_clock::now();

        FractionalScheme liftedScheme;

        int step = 0;
        bool reconstructed = schemes[i].reconstruct(liftedScheme) && liftedScheme.validate();

        if (!reconstructed) {
            auto lifter = schemes[i].toLift();

            while (step < steps && !reconstructed && lifter.lift()) {
                reconstructed = lifter.reconstruct(liftedScheme) && liftedScheme.validate();
                step++;
            }
        }

        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes[i] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0;
        std::string status;

        if (reconstructed) {
            if (canonize)
                liftedScheme.canonize();

            status = "reconstructed in " + liftedScheme.getRing();
            std::string path = getSavePath(liftedScheme, i, outputPath, format);

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

        std::stringstream ss;
        ss << "| " << std::setw(6) << (i + 1) << " | ";
        ss << std::setw(9) << schemes[i].getDimension() << " | ";
        ss << std::setw(4) << schemes[i].getRank() << " | ";
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
    int maxMatrixElements = std::stoi(parser["--int-width"]);

    if (maxMatrixElements <= 16)
        return runLiftSchemes<Scheme<uint16_t>>(parser);

    if (maxMatrixElements <= 32)
        return runLiftSchemes<Scheme<uint32_t>>(parser);

    if (maxMatrixElements <= 64)
        return runLiftSchemes<Scheme<uint64_t>>(parser);

    return runLiftSchemes<Scheme<__uint128_t>>(parser);
}

int main(int argc, char *argv[]) {
    ArgParser parser("lift", "Lift schemes from Z2/Z3 field to general");
    parser.addChoices("--ring", "-r", ArgType::String, "Coefficient ring: Z2 - {0, 1}, Z3 - {0, 1, 2}", {"Z2", "Z3"}, "", true);
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.addChoices("--format", "-f", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme(s)", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for lifted schemes", "schemes");
    parser.add("--multiple", "-m", ArgType::Flag, "Read multiple schemes from file, with total count on first line");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");

    parser.addSection("Lifting parameters");
    parser.add("--steps", "-k", ArgType::Natural, "Number of Hensel lifting steps", "10");
    parser.add("--canonize", "-c", ArgType::Flag, "Canonize reconstructed schemes");

    parser.addSection("Other parameters");
    parser.addChoices("--int-width", ArgType::String, "Integer bit width (16/32/64/128), determines maximum matrix elements", {"16", "32", "64", "128"}, "64");

    if (!parser.parse(argc, argv))
        return 0;

    if (parser["--ring"] == "Z2")
        return runLiftSchemesSizes<BinaryScheme>(parser);

    if (parser["--ring"] == "Z3")
        return runLiftSchemesSizes<Mod3Scheme>(parser);

    return 0;
}
