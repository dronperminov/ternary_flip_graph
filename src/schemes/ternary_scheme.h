#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cassert>
#include <algorithm>

#include "../entities/ternary_vector.hpp"
#include "../entities/flip_set.h"

typedef uint64_t vec_type;

class TernaryScheme {
    int dimension[3];
    int elements[3];
    int rank;
    std::vector<TernaryVector<vec_type>> uvw[3];
    FlipSet flips[3];
    std::vector<int> indices;

    std::uniform_int_distribution<int> boolDistribution;
    std::uniform_int_distribution<int> ijkDistribution;
public:
    TernaryScheme();
    TernaryScheme(const TernaryScheme &scheme);

    bool initializeNaive(int n1, int n2, int n3);
    bool read(const std::string &path);

    int getRank() const;
    int getComplexity() const;
    int getDimension(int index) const;
    int getAvailableFlips() const;

    bool tryFlip(std::mt19937 &generator);
    bool tryPlus(std::mt19937 &generator);
    bool trySplit(std::mt19937 &generator);
    bool tryExpand(std::mt19937 &generator);
    bool tryReduce();

    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void copy(const TernaryScheme &scheme);

    bool validate() const;
private:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const TernaryVector<vec_type> &u, const TernaryVector<vec_type> &v, const TernaryVector<vec_type> &w);

    void flip(int i, int j, int k, int index1, int index2);
    bool plus(int i, int j, int k, int index1, int index2, int variant);
    void split(int i, int j, int k, int index1, int index2);
    void reduceAdd(int i, int index1, int index2);
    void reduceSub(int i, int index1, int index2);
    bool checkFlipReduce(int j, int k, int index1, int index2);

    bool fixSigns();
    bool validateDimensions() const;
    bool validateEquation(int i, int j, int k) const;
    void saveMatrix(std::ofstream &f, std::string name, const std::vector<TernaryVector<vec_type>> &vectors) const;
};
