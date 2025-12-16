#include "flip_graph.h"

FlipGraph::FlipGraph(int count, const std::string outputPath, int threads, size_t flipIterations, size_t plusIterations, double reduceProbability, int seed, int topCount) : uniform(0.0, 1.0) {
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, count);
    this->flipIterations = flipIterations;
    this->plusIterations = plusIterations;
    this->resetIterations = 100000000;
    this->reduceProbability = reduceProbability;
    this->seed = seed;
    this->topCount = std::min(topCount, count);

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

    schemes.resize(count);
    schemesBest.resize(count);
    bestRanks.resize(count);
    flips.resize(count);
    indices.resize(count);
}

void FlipGraph::run(const TernaryScheme &scheme, int targetRank) {
    initialize(scheme);

    auto startTime = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    for (size_t iteration = 0; bestRank > targetRank; iteration++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        runIteration();
        updateBest(iteration);
        auto t2 = std::chrono::high_resolution_clock::now();
        elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);

        report(iteration + 1, startTime, elapsedTimes);
    }
}

void FlipGraph::initialize(const TernaryScheme &scheme) {
    bestRank = scheme.getRank();

    for (int i = 0; i < count; i++) {
        schemes[i].copy(scheme);
        schemesBest[i].copy(scheme);
        bestRanks[i] = bestRank;
        flips[i] = 0;
        indices[i] = i;
    }
}

void FlipGraph::runIteration()  {
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        randomWalk(schemes[i], schemesBest[i], flips[i], bestRanks[i], generators[omp_get_thread_num()]);
}

void FlipGraph::updateBest(size_t iteration) {
    std::partial_sort(indices.begin(), indices.begin() + topCount, indices.end(), [this](int index1, int index2) {
        return compare(index1, index2);
    });

    int top = indices[0];
    if (bestRanks[top] >= bestRank)
        return;

    if (!schemesBest[top].validate()) {
        std::cout << "Unable to save: scheme invalid" << std::endl;
        return;
    }

    std::string path = getSavePath(schemesBest[top], iteration, outputPath);
    schemesBest[top].saveJson(path + ".json");
    schemesBest[top].saveTxt(path + ".txt");

    std::cout << "Rank was improved from " << bestRank << " to " << bestRanks[top] << ", scheme was saved to \"" << path << "\"" << std::endl;
    bestRank = bestRanks[top];

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        if (uniform(generators[omp_get_thread_num()]) < 0.5)
            schemes[i].copy(schemesBest[indices[0]]);
}

void FlipGraph::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    std::cout << std::left;
    std::cout << "+-------------------------------------------------------------------+" << std::endl;
    std::cout << "| ";
    std::cout << "dimension: " << std::setw(8) << prettyDimension(schemes[indices[0]]) << "   ";
    std::cout << "flip iters: " << std::setw(7) << prettyInt(flipIterations) << "   ";
    std::cout << std::right << std::setw(21) << ("best rank: " + std::to_string(bestRank));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "threads: " << std::setw(10) << threads << "   ";
    std::cout << "plus iters: " << std::setw(7) << prettyInt(plusIterations) << "   ";
    std::cout << std::right << std::setw(21) << ("iteration: " + std::to_string(iteration));
    std::cout << " |" << std::endl;

    std::cout << "| " << std::left;
    std::cout << "count: " << std::setw(12) << count << "   ";
    std::cout << "seed: " << std::setw(13) << seed << "   ";
    std::cout << std::right << std::setw(21) << ("elapsed: " + prettyTime(elapsed));
    std::cout << " |" << std::endl;

    std::cout << "+===================================================================+" << std::endl;
    std::cout << "| run id | best | curr | complexity | flips count | flips available |" << std::endl;
    std::cout << "+--------+------+------+------------+-------------+-----------------+" << std::endl;
    std::cout << std::right;

    for (int i = 0; i < topCount; i++) {
        std::cout << "| ";
        std::cout << std::setw(6) << indices[i] << " | ";
        std::cout << std::setw(4) << schemesBest[indices[i]].getRank() << " | ";
        std::cout << std::setw(4) << schemes[indices[i]].getRank() << " | ";
        std::cout << std::setw(10) << schemes[indices[i]].getComplexity() << " | ";
        std::cout << std::setw(11) << prettyInt(flips[indices[i]]) << " | ";
        std::cout << std::setw(15) << schemes[indices[i]].getAvailableFlips() << " |";
        std::cout << std::endl;
    }

    std::cout << "+--------+------+------+------------+-------------+-----------------+" << std::endl;
    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

void FlipGraph::randomWalk(TernaryScheme &scheme, TernaryScheme &schemeBest, size_t &flipsCount, int &bestRank, std::mt19937 &generator) {
    for (size_t iteration = 0; iteration < flipIterations; iteration++) {
        int prevRank = scheme.getRank();

        if (!scheme.tryFlip(generator)) {
            if (scheme.tryExpand(generator))
                flipsCount = 0;

            continue;
        }

        if (uniform(generator) < reduceProbability && scheme.tryReduce())
            flipsCount = 0;

        int rank = scheme.getRank();
        if (rank < prevRank)
            flipsCount = 0;

        flipsCount++;

        if (rank < bestRank) {
            schemeBest.copy(scheme);
            bestRank = rank;
        }

        if (flipsCount >= plusIterations && rank < bestRank + 1 && scheme.tryExpand(generator)) {
            flipsCount = 0;
        }
        else if (flipsCount >= resetIterations) {
            scheme.copy(schemeBest);
            flipsCount = 0;
        }
    }
}

bool FlipGraph::compare(int index1, int index2) const {
    int bestRank1 = schemesBest[index1].getRank();
    int bestRank2 = schemesBest[index2].getRank();

    if (bestRank1 != bestRank2)
        return bestRank1 < bestRank2;

    int rank1 = schemes[index1].getRank();
    int rank2 = schemes[index2].getRank();

    if (rank1 != rank2)
        return rank1 < rank2;

    int complexity1 = schemes[index1].getComplexity();
    int complexity2 = schemes[index2].getComplexity();

    if (complexity1 != complexity2)
        return complexity1 < complexity2;

    return index1 < index2;
}

std::string FlipGraph::getSavePath(const TernaryScheme &scheme, int iteration, const std::string path) const {
    std::stringstream ss;
    ss << path << "/";
    ss << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2);
    ss << "_m" << scheme.getRank();
    ss << "_c" << scheme.getComplexity();
    ss << "_iteration" << iteration;
    return ss.str();
}

std::string FlipGraph::prettyInt(size_t value) const {
    std::stringstream ss;

    if (value < 1000)
        ss << value;
    else if (value < 1000000)
        ss << std::setprecision(2) << std::fixed << (value / 1000.0) << "K";
    else
        ss << std::setprecision(2) << std::fixed << (value / 1000000.0) << "M";

    return ss.str();
}

std::string FlipGraph::prettyDimension(const TernaryScheme &scheme) const {
    std::stringstream ss;
    ss << scheme.getDimension(0) << "x" << scheme.getDimension(1) << "x" << scheme.getDimension(2);
    return ss.str();
}

std::string FlipGraph::prettyTime(double elapsed) const {
    std::stringstream ss;

    if (elapsed < 60) {
        ss << std::setprecision(2) << std::fixed << elapsed;
    }
    else {
        int seconds = int(elapsed + 0.5);
        int hours = seconds / 3600;
        int minutes = (seconds % 3600) / 60;

        ss << std::setw(2) << std::setfill('0') << hours << ":";
        ss << std::setw(2) << std::setfill('0') << minutes << ":";
        ss << std::setw(2) << std::setfill('0') << (seconds % 60);
    }

    return ss.str();
}
