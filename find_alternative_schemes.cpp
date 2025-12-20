#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <random>
#include <ctime>
#include <unordered_set>

#include "src/entities/arg_parser.h"
#include "src/schemes/binary_scheme.hpp"
#include "src/schemes/ternary_scheme.hpp"

template <template<typename> typename Scheme, typename T>
std::string getSavePath(const Scheme<T> &scheme, const std::string &outputPath, int version) {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2);
    ss << "_m" << scheme.getRank();
    ss << "_v" << std::setw(6) << std::setfill('0') << version;
    ss << "_" << scheme.getRing();
    ss << ".json";
    return ss.str();
}

template <template<typename> typename Scheme, typename T>
int runFindAlternativeSchemes(const ArgParser &parser) {
    std::string inputPath = parser.get("-i");
    std::string outputPath = parser.get("-o");
    size_t maxCount = std::stoul(parser.get("--max-count"));
    int seed = std::stoi(parser.get("--seed"));

    if (seed == 0)
        seed = time(0);

    Scheme<T> scheme;
    if (!scheme.read(inputPath))
        return -1;

    std::cout << "Start finding alternative schemes" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << "- max count: " << maxCount << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- dimension: " << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2) << std::endl;
    std::cout << "- rank: " << scheme.getRank() << std::endl;
    std::cout << "- ring: " << scheme.getRing() << std::endl;
    std::cout << std::endl;

    std::mt19937 generator(seed);
    std::unordered_set<std::string> hashes;
    hashes.insert(scheme.getHash());

    size_t count = 0;

    std::cout << "+-----------+-------------+------------+" << std::endl;
    std::cout << "| iteration | alternative | complexity |" << std::endl;
    std::cout << "+-----------+-------------+------------+" << std::endl;

    for (size_t iteration = 1; count < maxCount && scheme.tryFlip(generator); iteration++) {
        std::string hash = scheme.getHash();
        if (hashes.find(hash) != hashes.end())
            continue;

        count++;
        std::cout << "| " << std::setw(9) << iteration << " | " << std::setw(11) << count << " | " << std::setw(10) << scheme.getComplexity() << " |" << std::endl;

        std::string path = getSavePath(scheme, outputPath, count);
        scheme.saveJson(path);
        hashes.insert(hash);
    }

    std::cout << "+-----------+-------------+------------+" << std::endl;
    return 0;
}

int main(int argc, char **argv) {
    ArgParser parser("find_alternative_schemes", "Find alternative fast matrix multiplication schemes using flip graph");

    parser.add("-i", ArgType::String, "PATH", "path to input file with initial scheme", "NULL");
    parser.add("-o", ArgType::String, "PATH", "output directory for alternative schemes", "schemes");
    parser.add("--ring", ArgType::String, "Z2/ZT", "coefficient ring: Z2 ({0, 1}) or ZT ({-1, 0, 1})", "ZT");
    parser.add("--max-count", ArgType::Natural, "INT", "number of alternative schemes", "10000");
    parser.add("--seed", ArgType::Natural, "INT", "random seed (0 = time-based)", "0");

    if (!parser.parse(argc, argv))
        return 0;

    std::string ring = parser.get("--ring");
    if (ring == "Z2" || ring == "binary")
        return runFindAlternativeSchemes<BinaryScheme, uint64_t>(parser);

    if (ring == "ZT" || ring == "ternary")
        return runFindAlternativeSchemes<TernaryScheme, uint64_t>(parser);

    return 0;
}
