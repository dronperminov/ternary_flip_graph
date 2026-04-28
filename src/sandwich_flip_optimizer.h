#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <omp.h>

#include "parameters/sandwiching_parameters.h"
#include "parameters/sandwich_flip_parameters.h"
#include "parameters/scale_parameters.h"
#include "parameters/plus_parameters.h"
#include "schemes/fractional_scheme.h"
#include "utils.h"

struct Weight {
    double omega;
    int flips;
    int flips3[3];
    int fractions;
    int fractions3[3];
    int denominator;
    int denominatorCount;
    int numerator;
    int numeratorCount;
    int weight;
    int complexity;
    double norm;
};


class SandwichFlipOptimizer {
    int count;
    std::string outputPath;
    int threads;

    SandwichFlipParameters sandwichFlipParameters;
    SandwichingParameters sandwichingParameters;
    ScaleParameters scaleParameters;
    PlusParameters plusParameters;

    int seed;
    bool noVerify;
    size_t maxImprovements;
    size_t improvementsIndex;
    std::string format;

    std::vector<FractionalScheme> schemes;
    std::vector<FractionalScheme> bestSchemes;
    std::vector<FractionalScheme> improvements;
    std::vector<Weight> weights;
    std::vector<Matrix> uvw[3];
    std::vector<Matrix> uvw1[3];
    Weight bestWeight;

    std::vector<std::mt19937> generators;
    std::uniform_real_distribution<double> uniform;
public:
    SandwichFlipOptimizer(int count, const std::string &outputPath, int threads, const SandwichFlipParameters &sandwichFlipParameters, const SandwichingParameters &sandwichingParameters, const ScaleParameters &scaleParameters, const PlusParameters &plusParameters, int seed, size_t maxImprovements, const std::string &format);

    bool initializeFromFile(const std::string &path, bool checkCorrectness, bool integer);

    void run(size_t maxNoImprovements);
private:
    void initialize();
    void optimize(FractionalScheme &scheme, FractionalScheme &bestScheme, Weight &weight, Matrix &u, Matrix &v, Matrix &w, Matrix &u1, Matrix &v1, Matrix &w1, std::mt19937 &generator);

    void resetImprovements();
    void addImprovement(const FractionalScheme &scheme);

    void makeSandwiching(FractionalScheme &scheme, Matrix &u, Matrix &v, Matrix &w, Matrix &u1, Matrix &v1, Matrix &w1, std::mt19937 &generator);
    void makeScale(FractionalScheme &scheme, std::mt19937 &generator);
    void makePlus(FractionalScheme &scheme, std::mt19937 &generator);

    void randomMatrix(Matrix &matrix, Matrix &inverse, int n, std::mt19937 &generator);
    void saveScheme(const FractionalScheme &scheme, const Weight &weight) const;
    void printHeader() const;

    Weight getWeight(const FractionalScheme &scheme, std::mt19937 &generator);
    bool compareWeight(const Weight &w1, const Weight &w2);
};
