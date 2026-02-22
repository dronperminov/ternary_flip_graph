#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <omp.h>


template <typename Scheme>
class SchemeOptimizer {
    int initialCount;
    int count;
    std::string outputPath;
    int threads;
    size_t flipIterations;
    double plusProbability;
    int plusDiff;
    int seed;
    double copyBestProbability;
    int goal;
    int topCount;
    std::string format;
    std::string metric;

    std::vector<Scheme> schemes;
    std::vector<Scheme> schemesBest;
    std::vector<int> bestMetrics;
    std::vector<int> indices;
    int bestMetric;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
public:
    SchemeOptimizer(int count, const std::string &outputPath, int threads, size_t flipIterations, double plusProbability, int pllusDiff, int seed, double copyBestProbability, const std::string &metric, bool maximize, int topCount, const std::string &format);

    bool initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness);
    void run(int maxNoImprovements);
private:
    void initialize();
    void optimizeIteration();

    void optimize(Scheme &scheme, Scheme &schemeBest, int &bestMetric, std::mt19937 &generator);
    bool updateBest();
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;

    int getMetric(const Scheme &scheme) const;
    std::string getSavePath(const Scheme &scheme) const;
};

template <typename Scheme>
SchemeOptimizer<Scheme>::SchemeOptimizer(int count, const std::string &outputPath, int threads, size_t flipIterations, double plusProbability, int plusDiff, int seed, double copyBestProbability, const std::string &metric, bool maximize, int topCount, const std::string &format) : uniform(0.0, 1.0) {
    this->initialCount = 0;
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, count);
    this->flipIterations = flipIterations;
    this->plusProbability = plusProbability;
    this->plusDiff = plusDiff;
    this->seed = seed;
    this->copyBestProbability = copyBestProbability;
    this->metric = metric;
    this->goal = maximize ? -1 : 1;
    this->topCount = std::min(topCount, count);
    this->format = format;
    this->bestMetric = 0;

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

    schemes.resize(count);
    schemesBest.resize(count);
    indices.resize(count);
    bestMetrics.resize(count);

    for (int i = 0; i < count; i++)
        indices[i] = i;
}

template <typename Scheme>
bool SchemeOptimizer<Scheme>::initializeFromFile(const std::string &path, bool multiple, bool checkCorrectness) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    initialCount = 1;

    if (multiple)
        f >> initialCount;

    std::cout << "Start reading " << std::min(initialCount, count) << " / " << initialCount << " schemes from \"" << path << "\"" << std::endl;

    bool valid = true;
    for (int i = 0; i < initialCount && i < count && valid; i++)
        valid = schemes[i].read(f, checkCorrectness);

    f.close();

    if (!valid)
        return false;

    #pragma omp parallel for num_threads(threads)
    for (int i = initialCount; i < count; i++)
        schemes[i].copy(schemes[i % initialCount]);

    return initialCount > 0;
}

template <typename Scheme>
void SchemeOptimizer<Scheme>::run(int maxNoImprovements) {
    initialize();

    auto startTime = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    int noImprovements = 0;

    for (size_t iteration = 1; noImprovements < maxNoImprovements; iteration++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        optimizeIteration();
        bool improved = updateBest();
        auto t2 = std::chrono::high_resolution_clock::now();

        elapsedTimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0);
        report(iteration, startTime, elapsedTimes);

        if (improved) {
            noImprovements = 0;
        }
        else {
            noImprovements++;
            std::cout << "No improvements for " << noImprovements << " iterations" << std::endl;
        }
    }
}

template <typename Scheme>
void SchemeOptimizer<Scheme>::initialize() {
    bestMetric = getMetric(schemes[0]);

    for (int i = 1; i < count && i < initialCount; i++)
        bestMetric = std::min(bestMetric, getMetric(schemes[i]));

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++) {
        bestMetrics[i] = bestMetric;
        schemesBest[i].copy(schemes[i]);
    }

    std::cout << "Initialized. Initial best " << metric << ": " << bestMetric << std::endl;
}

template <typename Scheme>
void SchemeOptimizer<Scheme>::optimizeIteration() {
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        optimize(schemes[i], schemesBest[i], bestMetrics[i], generators[omp_get_thread_num()]);
}

template <typename Scheme>
void SchemeOptimizer<Scheme>::optimize(Scheme &scheme, Scheme &schemeBest, int &bestMetric, std::mt19937 &generator) {
    int targetRank = schemeBest.getRank();

    for (size_t iteration = 0; iteration < flipIterations; iteration++) {
        if (!scheme.tryFlip(generator) || (scheme.getRank() < targetRank + plusDiff && uniform(generator) < plusProbability))
            scheme.tryExpand(generator);

        if (scheme.getRank() != targetRank)
            continue;

        int currMetric = getMetric(scheme);
        if ((currMetric - bestMetric) * goal < 0) {
            bestMetric = currMetric;
            schemeBest.copy(scheme);
        }
    }

    if (scheme.getRank() != targetRank)
        scheme.copy(schemeBest);
}

template <typename Scheme>
bool SchemeOptimizer<Scheme>::updateBest() {
    std::partial_sort(indices.begin(), indices.begin() + topCount, indices.end(), [this](int index1, int index2) {
        return (bestMetrics[index1] - bestMetrics[index2]) * goal < 0;
    });

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        if (uniform(generators[omp_get_thread_num()]) < copyBestProbability)
            schemes[i].copy(schemesBest[indices[0]]);

    int top = indices[0];
    if ((bestMetrics[top] - bestMetric) * goal >= 0)
        return false;

    if (!schemesBest[top].validate()) {
        std::cout << "Unable to save: scheme invalid" << std::endl;
        return false;
    }

    std::string path = getSavePath(schemesBest[top]);

    if (format == "json")
        schemesBest[top].saveJson(path);
    else
        schemesBest[top].saveTxt(path);

    std::cout << metric << " was improved from " << bestMetric << " to " << bestMetrics[top] << ", scheme was saved to \"" << path << "\"" << std::endl;
    bestMetric = bestMetrics[top];

    return true;
}

template <typename Scheme>
void SchemeOptimizer<Scheme>::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

    std::string metricHeader = "scheme " + metric;
    std::string left = std::string((23 - metricHeader.length()) / 2, ' ');
    std::string right = std::string((23 - metricHeader.length() - left.length()), ' ');

    std::cout << std::right;
    std::cout << "+----------------------------------+" << std::endl;
    std::cout << "| dimension       rank        ring |" << std::endl;
    std::cout << "| ";
    std::cout << std::setw(9) << schemes[indices[0]].getDimension() << "       ";
    std::cout << std::setw(4) << schemes[indices[0]].getRank() << "        ";
    std::cout << std::setw(4) << schemes[indices[0]].getRing();
    std::cout << " |" << std::endl;
    std::cout << "+----------------------------------+" << std::endl;
    std::cout << std::left;
    std::cout << "| count: " << std::setw(25) << (std::to_string(count) + " (" + std::to_string(threads) + " threads)") << " |" << std::endl;
    std::cout << "| seed: " << std::setw(26) << seed << " |" << std::endl;
    std::cout << "| best " << metric << ": " << std::setw(25 - metric.length()) << bestMetric << " |" << std::endl;
    std::cout << "| iteration: " << std::setw(21) << iteration << " |" << std::endl;
    std::cout << "| elapsed: " << std::setw(23) << prettyTime(elapsed) << " |" << std::endl;
    std::cout << "+==================================+" << std::endl;
    std::cout << "| runner | " << left << metricHeader << right << " |" << std::endl;
    std::cout << "|   id   |    best    |    curr    |" << std::endl;
    std::cout << "+--------+------------+------------+" << std::endl;

    for (int i = 0; i < topCount; i++) {
        int runner = indices[i];
        std::cout << "| ";
        std::cout << std::setw(6) << runner << " | ";
        std::cout << std::setw(10) << bestMetrics[runner] << " | ";
        std::cout << std::setw(10) << getMetric(schemes[runner]) << " | ";
        std::cout << std::endl;
    }

    std::cout << "+--------+------------+------------+" << std::endl;
    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
int SchemeOptimizer<Scheme>::getMetric(const Scheme &scheme) const {
    if (metric == "flips")
        return scheme.getAvailableFlips();

    return scheme.getComplexity();
}

template <typename Scheme>
std::string SchemeOptimizer<Scheme>::getSavePath(const Scheme &scheme) const {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_" << metric[0] << getMetric(scheme);
    ss << "_" << scheme.getRing();
    ss << "." << format;
    return ss.str();
}
