#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "src/utils.h"
#include "src/entities/arg_parser.h"
#include "src/schemes/fractional_scheme.h"

struct Weight {
    int flips;
    int fractions;
    int maxValue;
    int maxCount;
    int weight;
};

bool compareWeight(const Weight &w1, const Weight &w2, bool checkFlips) {
    if (checkFlips && w1.flips != w2.flips)
        return w1.flips > w2.flips;

    if (w1.fractions != w2.fractions)
        return w1.fractions < w2.fractions;

    if (w1.maxValue != w2.maxValue)
        return w1.maxValue < w2.maxValue;

    if (w1.maxCount != w2.maxCount)
        return w1.maxCount < w2.maxCount;
    
    return w1.weight < w2.weight;
}

Weight getWeight(const FractionalScheme &scheme) {
    Weight weight;
    weight.flips = scheme.getAvailableFlips();
    weight.fractions = scheme.getFractionsCount();
    weight.maxValue = scheme.getMaxAbsInteger();
    weight.maxCount = scheme.getAbsIntCount(weight.maxValue);
    weight.weight = scheme.getWeight();

    return weight;
}

std::ostream& operator<<(std::ostream &os, const Weight &weight) {
    os << "(";
    os << weight.flips << ", ";
    os << weight.fractions << ", ";
    os << weight.maxValue << ", ";
    os << weight.maxCount << ", ";
    os << weight.weight;
    os << ")";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> values) {
    os << "[";

    for (size_t i = 0; i < values.size(); i++)
        os << (i > 0 ? ", " : "") << values[i];

    os << "]";
    return os;
}

std::vector<Fraction> parseScales(const std::string &denominators) {
    std::vector<Fraction> scales = {Fraction(1)};

    std::stringstream ss(denominators);
    int denominator;

    while (ss >> denominator) {
        scales.push_back(Fraction(denominator, 1));
        scales.push_back(Fraction(1, denominator));
    }

    return scales;
}

void randomMatrix(Matrix &matrix, Matrix &inverse, int n, std::mt19937 &generator, int maxDen) {
    std::vector<int> values = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1};
    bool identity = generator() % 2 == 0;

    do {
        if (identity) {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    matrix[i * n + j] = i == j ? 1 : 0;
        }
        else {
            for (int i = 0; i < n * n; i++)
                matrix[i] = values[generator() % values.size()];
        }
    } while (!matrix.invertible(inverse) || inverse.maxDenominator() > maxDen);
}

void makeSandwiching(FractionalScheme &scheme, std::mt19937 &generator, int maxDen) {
    int n1 = scheme.getDimension(0);
    int n2 = scheme.getDimension(1);
    int n3 = scheme.getDimension(2);

    Matrix u(n1, n1), u1(n1, n1);
    Matrix v(n2, n2), v1(n2, n2);
    Matrix w(n3, n3), w1(n3, n3);

    randomMatrix(u, u1, n1, generator, maxDen);
    randomMatrix(v, v1, n2, generator, maxDen);
    randomMatrix(w, w1, n3, generator, maxDen);

    scheme.sandwiching(u, v, w, u1, v1, w1);
}

void makeScale(FractionalScheme &scheme, std::mt19937 &generator, const std::vector<Fraction> scales) {
    int index = generator() % scheme.getRank();
    Fraction alpha = scales[generator() % scales.size()];
    Fraction beta = scales[generator() % scales.size()];
    Fraction gamma = 1 / (alpha * beta);

    scheme.scale(index, alpha, beta, gamma);
}

void saveScheme(const FractionalScheme &scheme, const std::string &outputPath, const std::string &format) {
    std::stringstream ss;
    Weight weight = getWeight(scheme);
    ss << outputPath;
    ss << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_fr" << weight.fractions;
    ss << "_mv" << weight.maxValue;
    ss << "_mc" << weight.maxCount;
    ss << "_w" << weight.weight;
    ss << "_f" << weight.flips;
    ss << "_" << scheme.getRing();
    ss << "." << format;
    std::string path = ss.str();
    scheme.save(path);
}

void addImprovement(std::vector<FractionalScheme> &schemes, const FractionalScheme &scheme, size_t maxCount, size_t &position) {
    if (schemes.size() < maxCount)
        schemes.emplace_back(FractionalScheme());

    schemes[position].copy(scheme);
    position = (position + 1) % maxCount;
}

bool runSandwichFlip(const ArgParser &parser) {
    std::string inputPath = parser["--input-path"];
    std::string outputPath = parser["--output-path"];

    int minSteps = std::stoi(parser["--min-steps"]);
    int maxSteps = std::stoi(parser["--max-steps"]);
    int maxDenominator = std::stoi(parser["--max-denominator"]);

    size_t maxImprovements = parseNatural(parser["--max-improvements"]);
    size_t position = 0;

    double sandwichingProbability = std::stod(parser["--sandwiching-probability"]);
    double scaleProbability = std::stod(parser["--scale-probability"]);
    std::vector<Fraction> scales = parseScales(parser["--scale-denominators"]);

    int seed = std::stoi(parser["--seed"]);
    if (seed == 0)
        seed = time(0);

    std::string format = parser["--format"];
    bool verbose = parser.isSet("--verbose");
    bool fixFractions = parser.isSet("--fix-fractions");
    bool maximizeFlips = parser.isSet("--maximize-flips");

    std::cout << "Parsed parameters of the sandwich-flip tool:" << std::endl;
    std::cout << "- input path: " << inputPath << std::endl;
    std::cout << "- output path: " << outputPath << std::endl;
    std::cout << std::endl;
    std::cout << "Run parameters:" << std::endl;
    std::cout << "- sandwiching probability: " << sandwichingProbability << std::endl;
    std::cout << "- scale probability: " << scaleProbability << std::endl;
    if (scaleProbability > 0)
        std::cout << "- scale factors: " << scales << std::endl;
    std::cout << "- steps: " << minSteps << " .. " << maxSteps << std::endl;
    std::cout << "- max denominator: " << maxDenominator << std::endl;
    std::cout << "- max improvements: " << maxImprovements << std::endl;
    if (fixFractions)
        std::cout << "- fix fractions: yes" << std::endl;
    if (maximizeFlips)
        std::cout << "- maximize flips: yes" << std::endl;
    std::cout << std::endl;
    std::cout << "Other parameters:" << std::endl;
    std::cout << "- seed: " << seed << std::endl;
    std::cout << "- format: " << format << std::endl;
    std::cout << std::endl;

    if (!makeDirectory(outputPath))
        return false;

    double sumP = sandwichingProbability + scaleProbability;

    if (sumP > 1) {
        std::cout << "Invalid probabilities: sandwiching + scale > 1" << std::endl;
        return false;
    }

    FractionalScheme initScheme;
    if (!initScheme.read(inputPath, !parser.isSet("--no-verify"), parser.isSet("--integer")))
        return false;

    Weight bestWeight = getWeight(initScheme);
    std::cout << "Successfully read initial scheme:" << std::endl;
    std::cout << "- dimension: " << initScheme.getDimension() << std::endl;
    std::cout << "- rank: " << initScheme.getRank() << std::endl;
    std::cout << "- flips: " << initScheme.getAvailableFlips() << std::endl;
    std::cout << "- values: " << initScheme.getUniqueValues() << std::endl;
    std::cout << "- weight: " << bestWeight << std::endl;
    std::cout << std::endl;

    if (initScheme.getAvailableFlips() == 0 && sumP < 1) {
        sandwichingProbability /= sumP;
        scaleProbability /= sumP;

        std::cout << "Initial scheme has no flips, probabilities changed:" << std::endl;
        std::cout << "- sandwiching probability: " << sandwichingProbability << std::endl;
        std::cout << "- scale probability: " << scaleProbability << std::endl;
        std::cout << std::endl;
    }

    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    std::vector<FractionalScheme> improvements;
    addImprovement(improvements, initScheme, maxImprovements, position);

    FractionalScheme scheme;

    while (1) {
        int steps = minSteps + generator() % (maxSteps - minSteps + 1);
        scheme.copy(improvements[generator() % improvements.size()]);

        for (int step = 0; step < steps; step++) {
            double p = uniform(generator);

            if (p < sandwichingProbability) {
                makeSandwiching(scheme, generator, maxDenominator);
            }
            else if (p < scaleProbability) {
                makeScale(scheme, generator, scales);
            }
            else {
                scheme.tryFlip(generator);
            }

            if (fixFractions)
                scheme.fixFractions();

            Weight weight = getWeight(scheme);

            if (compareWeight(weight, bestWeight, maximizeFlips)) {
                if (!scheme.validate()) {
                    std::cout << "Invalid scheme, skip" << std::endl;
                    break;
                }

                bestWeight = weight;
                saveScheme(scheme, outputPath, format);
                addImprovement(improvements, scheme, maxImprovements, position);
                std::cout << "New improvement (step: " << (step + 1) << "): " << bestWeight << std::endl;
            }

            if (verbose)
                std::cout << weight << "    " << bestWeight << std::endl;
        }
    }

    return true;
}

int main(int argc, char *argv[]) {
    ArgParser parser("sandwich_flip", "Try to optimize scheme using random flips, scales and sandwiching");

    parser.addSection("Input / output");
    parser.add("--input-path", "-i", ArgType::Path, "Path to input file with initial scheme", "", true);
    parser.add("--output-path", "-o", ArgType::Path, "Output directory for improved schemes", "schemes");
    parser.add("--no-verify", ArgType::Flag, "Skip checking Brent equations for correctness");
    parser.add("--integer", ArgType::Flag, "Read scheme as integer");

    parser.addSection("Run parameters");
    parser.add("--min-steps", ArgType::Natural, "Minimum number of random steps", "1");
    parser.add("--max-steps", ArgType::Natural, "Maximum number of random steps", "20");
    parser.add("--max-denominator", ArgType::Natural, "Maximum denominator of sandwich matrices", "1");
    parser.add("--max-improvements", ArgType::Natural, "Number of last improved schemes", "8");
    parser.add("--sandwiching-probability", ArgType::Real, "Probability of sandwiching operation, from 0.0 to 1.0", "0.5");
    parser.add("--scale-probability", ArgType::Real, "Probability of scale operation, from 0.0 to 1.0", "0.0");
    parser.add("--scale-denominators", ArgType::String, "Denominators for scale operation (for example, \"2 3 5\")", "2");
    parser.add("--fix-fractions", ArgType::Flag, "Try to convert fractions to integers");
    parser.add("--maximize-flips", ArgType::Flag, "Check flips count during comparison first");

    parser.addSection("Other parameters");
    parser.add("--seed", ArgType::Natural, "Random seed, 0 uses time-based seed", "0");
    parser.addChoices("--format", ArgType::String, "Output format for saved schemes", {"json", "txt"}, "json");
    parser.add("--verbose", "-v", ArgType::Flag, "Show every step weight");

    if (!parser.parse(argc, argv))
        return 0;

    if (!runSandwichFlip(parser))
        return -1;

    return 0;
}
