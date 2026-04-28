#include "sandwich_flip_optimizer.h"

SandwichFlipOptimizer::SandwichFlipOptimizer(int count, const std::string &outputPath, int threads, const SandwichFlipParameters &sandwichFlipParameters, const SandwichingParameters &sandwichingParameters, const ScaleParameters &scaleParameters, const PlusParameters &plusParameters, int seed, size_t maxImprovements, const std::string &format) : uniform(0.0, 1.0) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = threads;

    this->sandwichFlipParameters = sandwichFlipParameters;
    this->sandwichingParameters = sandwichingParameters;
    this->scaleParameters = scaleParameters;
    this->plusParameters = plusParameters;

    this->seed = seed;
    this->maxImprovements = maxImprovements;
    this->format = format;

    generators = initRandomGenerators(seed, threads);
    schemes.resize(count);
    bestSchemes.resize(count);
    weights.resize(count);

    resetImprovements();
}

bool SandwichFlipOptimizer::initializeFromFile(const std::string &path, bool checkCorrectness, bool integer) {
    std::cout << "Start reading scheme from \"" << path << "\"" << std::endl;

    FractionalScheme scheme;
    if (!scheme.read(path, checkCorrectness, integer))
        return false;

    outputPath += "/" + scheme.getDimension();
    if (!makeDirectory(outputPath))
        return false;

    addImprovement(scheme);

    noVerify = !checkCorrectness;
    bestWeight = getWeight(scheme, generators[0]);

    std::cout << "Successfully read initial scheme:" << std::endl;
    std::cout << "- dimension: " << scheme.getDimension() << std::endl;
    std::cout << "- rank: " << scheme.getRank() << std::endl;
    std::cout << "- omega: " << std::setprecision(15) << std::fixed << scheme.getOmega() << std::endl;

    if (sandwichFlipParameters.minimizeNorm)
        std::cout << "- Frobenius norm: " << bestWeight.norm << std::endl;

    if (sandwichFlipParameters.minimizeOmega)
        std::cout << "- structure omega: " << std::setprecision(15) << std::fixed << bestWeight.omega << std::endl;

    std::cout << "- flips: " << bestWeight.flips << " (" << bestWeight.flips3[0] << " / " << bestWeight.flips3[1] << " / " << bestWeight.flips3[2] << ")" << std::endl;
    std::cout << "- values: " << scheme.getUniqueValues() << std::endl;
    std::cout << "- fractions: " << bestWeight.fractions << std::endl;
    std::cout << "- max denominator: " << bestWeight.denominator << " (count: " << bestWeight.denominatorCount << ")" << std::endl;
    std::cout << "- max abs value: " << bestWeight.numerator << " (count: " << bestWeight.numeratorCount << ")" << std::endl;
    std::cout << "- weight: " << bestWeight.weight << std::endl;
    std::cout << "- complexity: " << bestWeight.complexity << std::endl;
    std::cout << std::endl;

    double sumP = sandwichingParameters.probability + scaleParameters.probability;
    if (plusParameters.probability == 0 && scheme.getAvailableFlips() == 0 && sumP < 1) {
        sandwichingParameters.probability /= sumP;
        scaleParameters.probability /= sumP;

        std::cout << "Initial scheme has no flips, probabilities changed:" << std::endl;
        std::cout << "- sandwiching probability: " << sandwichingParameters.probability << std::endl;
        std::cout << "- scale probability: " << scaleParameters.probability << std::endl;
        std::cout << std::endl;
    }

    return true;
}

void SandwichFlipOptimizer::run(size_t maxNoImprovements) {
    initialize();
    printHeader();

    auto startTime = std::chrono::high_resolution_clock::now();

    for (size_t iteration = 0; iteration < maxNoImprovements; iteration++) {
        #pragma omp parallel for num_threads(threads) schedule(dynamic, 32)
        for (int i = 0; i < count; i++)
            optimize(schemes[i], bestSchemes[i], weights[i], uvw[0][i], uvw[1][i], uvw[2][i], uvw1[0][i], uvw1[1][i], uvw1[2][i], generators[omp_get_thread_num()]);

        int imin = -1;
        for (int i = 0; i < count; i++) {
            if (compareWeight(weights[i], bestWeight)) {
                bestWeight = weights[i];
                imin = i;
            }
        }

        if (imin > -1) {
            iteration = 0;

            addImprovement(bestSchemes[imin]);
            saveScheme(bestSchemes[imin], bestWeight);

            double elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

            std::cout << "| ";
            if (sandwichFlipParameters.minimizeNorm)
                std::cout << std::setw(20) << std::setprecision(15) << bestWeight.norm << " | ";

            if (sandwichFlipParameters.minimizeOmega)
                std::cout << std::setw(17) << std::setprecision(15) << std::fixed << bestWeight.omega << " | ";

            std::cout << std::setw(5) << bestWeight.flips << " | ";
            std::cout << std::setw(3) << bestWeight.flips3[0] << " | ";
            std::cout << std::setw(3) << bestWeight.flips3[1] << " | ";
            std::cout << std::setw(3) << bestWeight.flips3[2] << " | ";
            std::cout << std::setw(9) << bestWeight.fractions << " | ";
            std::cout << std::setw(7) << bestWeight.denominator << " | ";
            std::cout << std::setw(7) << bestWeight.denominatorCount << " | ";
            std::cout << std::setw(7) << bestWeight.numerator << " | ";
            std::cout << std::setw(7) << bestWeight.numeratorCount << " | ";
            std::cout << std::setw(10) << bestWeight.weight << " | ";
            std::cout << std::setw(10) << bestWeight.complexity << " | ";
            std::cout << std::setw(11) << prettyTime(elapsedTime) << " |";
            std::cout << std::endl;
        }
    }
}

void SandwichFlipOptimizer::initialize() {
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        weights[i] = bestWeight;

    int n1 = improvements[0].getDimension(0);
    int n2 = improvements[0].getDimension(1);
    int n3 = improvements[0].getDimension(2);

    for (int i = 0; i < count; i++) {
        uvw[0].push_back(Matrix(n1, n1));
        uvw1[0].push_back(Matrix(n1, n1));

        uvw[1].push_back(Matrix(n2, n2));
        uvw1[1].push_back(Matrix(n2, n2));

        uvw[2].push_back(Matrix(n3, n3));
        uvw1[2].push_back(Matrix(n3, n3));
    }
}

void SandwichFlipOptimizer::optimize(FractionalScheme &scheme, FractionalScheme &bestScheme, Weight &weight, Matrix &u, Matrix &v, Matrix &w, Matrix &u1, Matrix &v1, Matrix &w1, std::mt19937 &generator) {
    scheme.copy(improvements[generator() % improvements.size()]);

    int steps = sandwichFlipParameters.minSteps + generator() % (sandwichFlipParameters.maxSteps - sandwichFlipParameters.minSteps + 1);
    int rank = scheme.getRank();

    for (int step = 0; step < steps; step++) {
        double p = uniform(generator);

        if (p < sandwichingParameters.probability) {
            makeSandwiching(scheme, u, v, w, u1, v1, w1, generator);
        }
        else if (p < sandwichingParameters.probability + scaleParameters.probability) {
            makeScale(scheme, generator);
        }
        else {
            if (scheme.getRank() == rank && uniform(generator) < plusParameters.probability) {
                makePlus(scheme, generator);
            }
            else {
                scheme.tryFlip(generator);
            }
        }

        if (scheme.getRank() != rank)
            continue;

        if (sandwichFlipParameters.fixFractions)
            scheme.fixFractions();

        Weight currWeight = getWeight(scheme, generator);

        if (compareWeight(currWeight, weight)) {
            if (!noVerify && !scheme.validateParallel())
                break;

            weight = currWeight;
            bestScheme.copy(scheme);
        }
    }
}

void SandwichFlipOptimizer::resetImprovements() {
    improvementsIndex = 0;
    improvements.clear();
}

void SandwichFlipOptimizer::addImprovement(const FractionalScheme &scheme) {
    if (improvements.size() < maxImprovements)
        improvements.emplace_back(FractionalScheme());

    improvements[improvementsIndex].copy(scheme);
    improvementsIndex = (improvementsIndex + 1) % maxImprovements;
}

void SandwichFlipOptimizer::makeSandwiching(FractionalScheme &scheme, Matrix &u, Matrix &v, Matrix &w, Matrix &u1, Matrix &v1, Matrix &w1, std::mt19937 &generator) {
    int n1 = scheme.getDimension(0);
    int n2 = scheme.getDimension(1);
    int n3 = scheme.getDimension(2);
    int variant = generator() % 100;

    u.identity();
    u1.identity();

    v.identity();
    v1.identity();

    w.identity();
    w1.identity();

    if (variant < 24) {
        randomMatrix(u, u1, n1, generator);
    }
    else if (variant < 48) {
        randomMatrix(v, v1, n2, generator);
    }
    else if (variant < 72) {
        randomMatrix(w, w1, n3, generator);
    }
    else if (variant < 81) {
        randomMatrix(u, u1, n1, generator);
        randomMatrix(v, v1, n2, generator);
    }
    else if (variant < 90) {
        randomMatrix(u, u1, n1, generator);
        randomMatrix(w, w1, n3, generator);
    }
    else if (variant < 99) {
        randomMatrix(v, v1, n2, generator);
        randomMatrix(w, w1, n3, generator);
    }
    else {
        randomMatrix(u, u1, n1, generator);
        randomMatrix(v, v1, n2, generator);
        randomMatrix(w, w1, n3, generator);
    }

    scheme.sandwiching(u, v, w, u1, v1, w1);
}

void SandwichFlipOptimizer::makeScale(FractionalScheme &scheme, std::mt19937 &generator) {
    Fraction alpha = scaleParameters.values[generator() % scaleParameters.values.size()];
    Fraction beta = scaleParameters.values[generator() % scaleParameters.values.size()];
    Fraction gamma = scaleParameters.values[generator() % scaleParameters.values.size()];

    int variant = generator() % 3;

    if (variant == 0) {
        gamma = 1 / (alpha * beta);
    }
    else if (variant == 1) {
        beta = 1 / (alpha * gamma);
    }
    else {
        alpha = 1 / (beta * gamma);
    }

    if (uniform(generator) < scaleParameters.fullProbability) {
        for (int index = 0; index < scheme.getRank(); index++)
            scheme.scale(index, alpha, beta, gamma);
    }
    else if (scaleParameters.maxRows > 1) {
        int rank = scheme.getRank();
        int rowsCount = 1 + generator() % std::min(scaleParameters.maxRows, rank);

        std::vector<int> rows(rank);
        for (int i = 0; i < rank; i++)
            rows[i] = i;
        std::shuffle(rows.begin(), rows.end(), generator);

        for (int i = 0; i < rowsCount; i++)
            scheme.scale(rows[i], alpha, beta, gamma);
    }
    else {
        scheme.scale(generator() % scheme.getRank(), alpha, beta, gamma);
    }
}

void SandwichFlipOptimizer::makePlus(FractionalScheme &scheme, std::mt19937 &generator) {
    int rank = scheme.getRank();

    if (generator() % 2) {
        scheme.plus(generator);
    }
    else {
        scheme.split(generator, plusParameters.values);
    }

    for (int j = 0; j < plusParameters.iterations && scheme.getRank() != rank; j++)
        scheme.tryFlip(generator);
}

void SandwichFlipOptimizer::randomMatrix(Matrix &matrix, Matrix &inverse, int n, std::mt19937 &generator) {
    std::vector<int> rows(n);
    std::vector<int> columns(n);

    for (int i = 0; i < n; i++) {
        rows[i] = i;
        columns[i] = i;
    }

    do {
        int variant = generator() % 10;

        if (variant < 9) {
            for (int i = 0; i < n*n; i++)
                matrix[i] = 0;

            std::shuffle(rows.begin(), rows.end(), generator);
            int rowsCount = sandwichingParameters.minRows + generator() % (sandwichingParameters.maxRows - sandwichingParameters.minRows + 1);

            for (int i = 0; i < rowsCount && i < n; i++) {
                std::shuffle(columns.begin(), columns.end(), generator);
                int columnsCount = sandwichingParameters.minColumns + generator() % (sandwichingParameters.maxColumns - sandwichingParameters.minColumns + 1);

                for (int j = 0; j < columnsCount && j < n; j++)
                    matrix(rows[i], columns[j]) = sandwichingParameters.nonZero[generator() % sandwichingParameters.nonZero.size()];
            }

            for (int i = rowsCount; i < n; i++)
                matrix(rows[i], rows[i]) = 1;

            if (generator() % 2)
                matrix.transpose();
        }
        else {
            for (int i = 0; i < n * n; i++)
                matrix[i] = sandwichingParameters.values[generator() % sandwichingParameters.values.size()];
        }
    } while (!matrix.invertible(inverse) || inverse.maxDenominator() > sandwichingParameters.maxDenominator);
}

void SandwichFlipOptimizer::saveScheme(const FractionalScheme &scheme, const Weight &weight) const {
    std::stringstream ss;
    ss << outputPath;
    ss << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();

    if (sandwichFlipParameters.minimizeNorm)
        ss << "_n" << std::fixed << std::setprecision(4) << weight.norm;

    if (sandwichFlipParameters.minimizeOmega)
        ss << "_o" << std::fixed << std::setprecision(15) << weight.omega;

    if (sandwichFlipParameters.maximizeFlips)
        ss << "_f" << weight.flips << "_" << weight.flips3[0] << "_" << weight.flips3[1] << "_" << weight.flips3[2];

    ss << "_md" << weight.denominator;
    ss << "_dc" << weight.denominatorCount;
    ss << "_fr" << weight.fractions << "_" << weight.fractions3[0] << "_" << weight.fractions3[1] << "_" << weight.fractions3[2];;
    ss << "_mn" << weight.numerator;
    ss << "_nc" << weight.numeratorCount;
    ss << "_w" << weight.weight;
    ss << "_c" << weight.complexity;

    if (!sandwichFlipParameters.maximizeFlips)
        ss << "_f" << weight.flips << "_" << weight.flips3[0] << "_" << weight.flips3[1] << "_" << weight.flips3[2];

    ss << "_" << scheme.getRing();
    ss << "." << format;

    std::string path = ss.str();
    scheme.save(path);
}

void SandwichFlipOptimizer::printHeader() const {
    std::vector<std::string> headers1;
    std::vector<std::string> headers2;

    if (sandwichFlipParameters.minimizeNorm) {
        headers1.push_back("      Frobenius     ");
        headers2.push_back("        norm        ");
    }

    if (sandwichFlipParameters.minimizeOmega) {
        headers1.push_back("    structure    ");
        headers2.push_back("      omega      ");
    }

    headers1.push_back("    available flips    ");
    headers2.push_back("total |  u  |  v  |  w ");

    headers1.push_back("fractions");
    headers2.push_back("  count  ");

    headers1.push_back(" max denominator ");
    headers2.push_back(" value  |  count ");

    headers1.push_back("max abs numerator");
    headers2.push_back(" value  |  count ");

    headers1.push_back("  weight  ");
    headers2.push_back("          ");

    headers1.push_back("  naive   ");
    headers2.push_back("complexity");

    headers1.push_back("  elapsed  ");
    headers2.push_back("    time   ");

    std::string divider = "+";
    std::string header1 = "|";
    std::string header2 = "|";

    for (size_t i = 0; i < headers1.size(); i++) {
        divider += std::string(headers1[i].size() + 2, '-') + "+";
        header1 += " " + headers1[i] + " |";
        header2 += " " + headers2[i] + " |";
    }

    std::cout << divider << std::endl;
    std::cout << header1 << std::endl;
    std::cout << header2 << std::endl;
    std::cout << divider << std::endl;
}

Weight SandwichFlipOptimizer::getWeight(const FractionalScheme &scheme, std::mt19937 &generator) {
    Weight weight;
    weight.norm = sandwichFlipParameters.minimizeNorm ? scheme.getFrobeniusNorm() : 0;
    weight.omega = sandwichFlipParameters.minimizeOmega ? scheme.getOptimalStructure(generator, 100, 1e-15).omega : scheme.getOmega();

    weight.flips = 0;
    weight.fractions = 0;
    for (int i = 0; i < 3; i++) {
        weight.flips3[i] = scheme.getAvailableFlips(i);
        weight.flips += weight.flips3[i];

        weight.fractions3[i] = scheme.getFractionsCount(i);
        weight.fractions += weight.fractions3[i];
    }

    weight.denominator = scheme.getMaxDenominator();
    weight.denominatorCount = scheme.getDenominatorCount(weight.denominator);
    weight.numerator = scheme.getMaxAbsNumerator();
    weight.numeratorCount = scheme.getAbsNumeratorCount(weight.numerator);
    weight.weight = scheme.getWeight();
    weight.complexity = scheme.getComplexity();
    return weight;
}

bool SandwichFlipOptimizer::compareWeight(const Weight &w1, const Weight &w2) {
    if (w1.fractions > sandwichFlipParameters.maxFractions || w1.denominator > sandwichFlipParameters.maxDenominator || w1.numerator > sandwichFlipParameters.maxNumerator || w1.weight > sandwichFlipParameters.maxWeight)
        return false;

    if (sandwichFlipParameters.minimizeNorm && w1.norm != w2.norm)
        return w1.norm < w2.norm;

    if (sandwichFlipParameters.minimizeOmega && w1.omega != w2.omega)
        return w1.omega < w2.omega;

    if (sandwichFlipParameters.maximizeFlips && w1.flips != w2.flips)
        return w1.flips > w2.flips;

    for (const char& c : sandwichFlipParameters.check) {
        if (c == 'd') {
            if (w1.denominator != w2.denominator)
                return w1.denominator < w2.denominator;

            if (w1.denominatorCount != w2.denominatorCount)
                return w1.denominatorCount < w2.denominatorCount;
        }

        if (c == 'f') {
            if (w1.fractions != w2.fractions)
                return w1.fractions < w2.fractions;

            if (w2.fractions == 0)
                continue;

            int min1 = std::min(w1.fractions3[0], std::min(w1.fractions3[1], w1.fractions3[2]));
            int min2 = std::min(w2.fractions3[0], std::min(w2.fractions3[1], w2.fractions3[2]));
            if (min1 != min2)
                return min1 < min2;

            int count1 = 0;
            int count2 = 0;
            for (int i = 0; i < 3; i++) {
                count1 += w1.fractions3[i] == 0;
                count2 += w2.fractions3[i] == 0;
            }

            if (count1 != count2)
                return count1 > count2;
        }

        if (c == 'n') {
            if (w1.numerator != w2.numerator)
                return w1.numerator < w2.numerator;

            if (w1.numeratorCount != w2.numeratorCount)
                return w1.numeratorCount < w2.numeratorCount;
        }

        if (c == 'w' && w1.weight != w2.weight)
            return w1.weight < w2.weight;

        if (c == 'c' && w1.complexity != w2.complexity)
            return w1.complexity < w2.complexity;
    }

    return false;
}
