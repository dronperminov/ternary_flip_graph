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
class ComplexityMinimizer {
    int initialCount;
    int count;
    std::string outputPath;
    int threads;
    size_t flipIterations;
    double plusProbability;
    int seed;
    int topCount;
    std::string format;

    std::vector<Scheme> schemes;
    std::vector<Scheme> schemesBest;
    std::vector<int> bestComplexities;
    std::vector<int> indices;
    int bestComplexity;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
public:
    ComplexityMinimizer(int count, const std::string &outputPath, int threads, size_t flipIterations, double plusProbability, int seed, int topCount, const std::string &format);

    bool initializeFromFile(const std::string &path);
    void run(int maxNoImprovements);
private:
    void initialize();
    void minimizeIteration();

    void minimize(Scheme &scheme, Scheme &schemeBest, int &bestComplexity, std::mt19937 &generator);
    bool updateBest();
    void report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const;

    std::string getSavePath(const Scheme &scheme) const;
};

template <typename Scheme>
ComplexityMinimizer<Scheme>::ComplexityMinimizer(int count, const std::string &outputPath, int threads, size_t flipIterations, double plusProbability, int seed, int topCount, const std::string &format) : uniform(0.0, 1.0) {
    this->initialCount = 0;
    this->count = count;
    this->outputPath = outputPath;
    this->threads = std::min(threads, count);
    this->flipIterations = flipIterations;
    this->plusProbability = plusProbability;
    this->seed = seed;
    this->topCount = std::min(topCount, count);
    this->format = format;
    this->bestComplexity = 0;

    for (int i = 0; i < threads; i++)
        generators.emplace_back(seed + i);

    schemes.resize(count);
    schemesBest.resize(count);
    indices.resize(count);
    bestComplexities.resize(count);

    for (int i = 0; i < count; i++)
        indices[i] = i;
}

template <typename Scheme>
bool ComplexityMinimizer<Scheme>::initializeFromFile(const std::string &path) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    f >> initialCount;

    std::cout << "Start reading " << std::min(initialCount, count) << " / " << initialCount << " schemes from \"" << path << "\"" << std::endl;

    bool valid = true;
    for (int i = 0; i < initialCount && i < count && valid; i++)
        valid = schemes[i].read(f);

    f.close();

    if (!valid)
        return false;

    #pragma omp parallel for num_threads(threads)
    for (int i = initialCount; i < count; i++)
        schemes[i].copy(schemes[i % initialCount]);

    return initialCount > 0;
}

template <typename Scheme>
void ComplexityMinimizer<Scheme>::run(int maxNoImprovements) {
    initialize();

    auto startTime = std::chrono::high_resolution_clock::now();
    std::vector<double> elapsedTimes;

    int noImprovements = 0;

    for (size_t iteration = 1; noImprovements < maxNoImprovements; iteration++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        minimizeIteration();
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
void ComplexityMinimizer<Scheme>::initialize() {
    bestComplexity = schemes[0].getComplexity();

    for (int i = 1; i < count && i < initialCount; i++)
        bestComplexity = std::min(bestComplexity, schemes[i].getComplexity());

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++) {
        bestComplexities[i] = bestComplexity;
        schemesBest[i].copy(schemes[i]);
    }

    std::cout << "Initialized. Initial best complexity: " << bestComplexity << std::endl;
}

template <typename Scheme>
void ComplexityMinimizer<Scheme>::minimizeIteration() {
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        minimize(schemes[i], schemesBest[i], bestComplexities[i], generators[omp_get_thread_num()]);
}

template <typename Scheme>
void ComplexityMinimizer<Scheme>::minimize(Scheme &scheme, Scheme &schemeBest, int &bestComplexity, std::mt19937 &generator) {
    for (size_t iteration = 0; iteration < flipIterations; iteration++) {
        if (!scheme.tryFlip(generator))
            break;

        if (scheme.getRank() != schemeBest.getRank())
            continue;

        int complexity = scheme.getComplexity();
        if (complexity < bestComplexity) {
            bestComplexity = complexity;
            schemeBest.copy(scheme);
        }

        if (uniform(generator) < plusProbability)
            scheme.tryPlus(generator);
    }
}

template <typename Scheme>
bool ComplexityMinimizer<Scheme>::updateBest() {
    std::partial_sort(indices.begin(), indices.begin() + topCount, indices.end(), [this](int index1, int index2) {
        return bestComplexities[index1] < bestComplexities[index2];
    });

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < count; i++)
        if (uniform(generators[omp_get_thread_num()]) < 0.5)
            schemes[i].copy(schemesBest[indices[0]]);

    int top = indices[0];
    if (bestComplexities[top] >= bestComplexity)
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

    std::cout << "Naive complexity was improved from " << bestComplexity << " to " << bestComplexities[top] << ", scheme was saved to \"" << path << "\"" << std::endl;
    bestComplexity = bestComplexities[top];

    return true;
}

template <typename Scheme>
void ComplexityMinimizer<Scheme>::report(size_t iteration, std::chrono::high_resolution_clock::time_point startTime, const std::vector<double> &elapsedTimes) const {
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

    double lastTime = elapsedTimes[elapsedTimes.size() - 1];
    double minTime = *std::min_element(elapsedTimes.begin(), elapsedTimes.end());
    double maxTime = *std::max_element(elapsedTimes.begin(), elapsedTimes.end());
    double meanTime = std::accumulate(elapsedTimes.begin(), elapsedTimes.end(), 0.0) / elapsedTimes.size();

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
    std::cout << "| best complexity: " << std::setw(15) << bestComplexity << " |" << std::endl;
    std::cout << "| iteration: " << std::setw(21) << iteration << " |" << std::endl;
    std::cout << "| elapsed: " << std::setw(23) << prettyTime(elapsed) << " |" << std::endl;
    std::cout << "+==================================+" << std::endl;
    std::cout << "| runner | naive scheme complexity |" << std::endl;
    std::cout << "|   id   |    best    |    curr    |" << std::endl;
    std::cout << "+--------+------------+------------+" << std::endl;

    for (int i = 0; i < topCount; i++) {
        int runner = indices[i];
        std::cout << "| ";
        std::cout << std::setw(6) << runner << " | ";
        std::cout << std::setw(10) << bestComplexities[runner] << " | ";
        std::cout << std::setw(10) << schemes[runner].getComplexity() << " | ";
        std::cout << std::endl;
    }

    std::cout << "+--------+------------+------------+" << std::endl;
    std::cout << "- iteration time (last / min / max / mean): " << prettyTime(lastTime) << " / " << prettyTime(minTime) << " / " << prettyTime(maxTime) << " / " << prettyTime(meanTime) << std::endl;
    std::cout << std::endl;
}

template <typename Scheme>
std::string ComplexityMinimizer<Scheme>::getSavePath(const Scheme &scheme) const {
    std::stringstream ss;
    ss << outputPath << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_c" << scheme.getComplexity();
    ss << "_" << scheme.getRing();
    ss << "." << format;
    return ss.str();
}
