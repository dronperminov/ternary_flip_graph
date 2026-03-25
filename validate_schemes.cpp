#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "src/entities/arg_parser.h"
#include "src/schemes/fractional_scheme.h"
#include "src/utils.h"

void showSchemeParameters(const FractionalScheme &scheme, bool showRing, bool showCoefficients, double elapsed) {
    std::cout << " (" << scheme.getDimension(0) << ", " << scheme.getDimension(1) << ", " << scheme.getDimension(2) << ": " << scheme.getRank() << ")";

    if (showRing)
        std::cout << ", ring: " << scheme.getRing();

    if (showCoefficients)
        std::cout << ", values: " << scheme.getUniqueValues();

    std::cout << ", elapsed: " << prettyTime(elapsed) << std::endl;
}

void validateSchemes(const std::string &path, bool multiple, bool showRing, bool showCoefficients, bool integer) {
    std::ifstream f(path);
    if (!f) {
        std::cout << "Unable to open file \"" << path << "\"" << std::endl;
        return;
    }

    int count = 1;
    if (multiple)
        f >> count;

    std::cout << "Start checking " << count << " schemes in \"" << path << "\"" << std::endl;

    int invalid = 0;

    for (int i = 0; i < count; i++) {
        auto startTime = std::chrono::high_resolution_clock::now();

        FractionalScheme scheme;
        if (!scheme.read(f, false, integer) || !scheme.validateParallel()) {
            std::cout << "- invalid scheme " << (i + 1) << " / " << count;
            invalid++;
        }
        else {
            std::cout << "- correct scheme " << (i + 1) << " / " << count;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() / 1000.0;

        showSchemeParameters(scheme, showRing, showCoefficients, elapsed);
    }

    f.close();

    std::cout << std::endl;
    if (count == 1) {
        std::cout << "Readed scheme is " << (invalid ? "invalid" : "correct") << std::endl;
    }
    else if (invalid == 0) {
        std::cout << "All " << count << " schemes are correct" << std::endl;
    }
    else {
        std::cout << invalid << " of " << count << " schemes are invalid" << std::endl;
    }
}

int main(int argc, char **argv) {
    ArgParser parser("validate_schemes", "Check validity of scheme(s)");
    parser.add("--input-path", "-i", ArgType::String, "Path to file with scheme(s)", "", true);
    parser.add("--multiple", "-m", ArgType::Flag, "Read multiple schemes from file, with total count on first line");
    parser.add("--show-ring", "-sr", ArgType::Flag, "Show the coefficient ring of checked schemes");
    parser.add("--show-coefficients", "-sc", ArgType::Flag, "Show the coefficient set of checked schemes");
    parser.addChoices("--format", "-f", ArgType::String, "Input scheme format", {"int", "frac"}, "frac", true);

    if (!parser.parse(argc, argv))
        return -1;

    std::string path = parser["--input-path"];
    bool multiple = parser.isSet("--multiple");
    bool showRing = parser.isSet("--show-ring");
    bool showCoefficients = parser.isSet("--show-coefficients");
    bool integer = parser["--format"] == "int";

    validateSchemes(path, multiple, showRing, showCoefficients, integer);
    return 0;
}
