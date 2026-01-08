#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/schemes/binary_scheme.hpp"
#include "src/schemes/fractional_scheme.hpp"

bool readSchemes(const std::string &inputPath, std::vector<BinaryScheme<uint64_t>> &schemes, bool multiple) {
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
        valid = schemes[i].read(f);

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

int main(int argc, char *argv[]) {
    ArgParser parser("lift", "Lift schemes from Z2/Z3 field to general");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme(s)", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for discovered schemes", "schemes");
    parser.add("--multiple", "-m", ArgType::Flag, "Read multiple schemes from file, with total count on first line");

    parser.addSection("Lifting parameters");
    parser.add("--steps", "-k", ArgType::Natural, "Number of Hensel lifting steps", "10");
    parser.add("--bound", "-b", ArgType::Natural, "Bound of rational reconstruction");
    parser.add("--canonize", "-c", ArgType::Flag, "Canonize reconstructed schemes");

    parser.addSection("Run parameters");
    parser.add("--threads", "-t", ArgType::Natural, "Number of OpenMP threads", std::to_string(omp_get_max_threads()));
    parser.addChoices("--format", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");

    if (!parser.parse(argc, argv))
        return 0;

    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    int steps = std::stoi(parser["--steps"]);
    int64_t mod = int64_t(1) << (steps + 1);
    int64_t bound = parser.isSet("--bound") ? std::stoll(parser["--bound"]) : (int64_t) std::sqrt(mod / 2.0);
    bool canonize = parser.isSet("--canonize");

    int threads = std::stoi(parser["--threads"]);
    std::string format = parser["--format"];

    if (!makeDirectory(outputPath))
        return -1;

    std::vector<BinaryScheme<uint64_t>> schemes;
    if (!readSchemes(inputPath, schemes, parser.isSet("--multiple")))
        return -1;

    std::cout << "Successfully read " << schemes.size() << " schemes from \"" << inputPath << "\"" << std::endl;
    std::cout << "Start " << steps << " steps Hensel lifting" << std::endl;
    std::cout << "- mod: " << mod << std::endl;
    std::cout << "- bound: " << bound << std::endl;
    std::cout << std::endl;

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < schemes.size(); i++) {
        int n1 = schemes[i].getDimension(0);
        int n2 = schemes[i].getDimension(1);
        int n3 = schemes[i].getDimension(2);
        int rank = schemes[i].getRank();

        std::vector<uint64_t> u;
        std::vector<uint64_t> v;
        std::vector<uint64_t> w;

        if (!schemes[i].lift(steps, u, v, w)) {
            std::cout << (i + 1) << ". Unable to lift scheme" << std::endl;
            continue;
        }

        FractionalScheme lifted;
        if (!lifted.reconstruct(n1, n2, n3, rank, u, v, w, mod, bound)) {
            std::cout << (i + 1) << ". Unable to make rational reconstruction" << std::endl;
            continue;
        }

        if (canonize)
            lifted.canonize();

        if (!lifted.validate()) {
            std::cout << (i + 1) << ". Reconstructed scheme is not valid" << std::endl;
            continue;
        }

        std::cout << (i + 1) << ". Successfully reconstructed scheme in " << lifted.getRing() << std::endl;
        std::string path = getSavePath(lifted, i, outputPath, format);

        if (format == "txt")
            lifted.saveTxt(path);
        else
            lifted.saveJson(path);
    }

    return 0;
}
