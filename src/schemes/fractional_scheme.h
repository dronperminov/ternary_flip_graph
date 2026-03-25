#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include "../algebra/fraction.h"
#include "../algebra/matrix.h"
#include "../entities/ranks.h"
#include "../entities/flip_structure_optimizer.h"
#include "base_scheme.h"

class FractionalScheme : public BaseScheme {
protected:
    std::vector<Fraction> uvw[3];
    FlipSet flipsNeg[3];
public:
    bool reconstruct(int n1, int n2, int n3, int rank, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, int64_t mod, int64_t bound);
    bool validate() const;
    bool validateParallel() const;

    bool read(const std::string &path, bool checkCorrectness, bool integer);
    bool read(std::istream &is, bool checkCorrectness, bool integer);

    bool isInteger() const;
    bool isTernary() const;

    int getAvailableFlips() const;
    int getFractionsCount() const;
    int getComplexity() const;
    int64_t getWeight() const;
    int getMaxAbsInteger() const;
    int getAbsIntCount(int value) const;
    int getMaxDenominator() const;
    std::string getRing() const;
    std::string getUniqueValues() const;

    FlipStructure getOptimalStructure(std::mt19937 &generator, int iterations = 250, double eps = 1e-15) const;

    std::string getTypeInvariant() const;

    bool tryFlip(std::mt19937 &generator);
    void plus(std::mt19937 &generator);
    void split(std::mt19937 &generator, const std::vector<Fraction> &values);
    void sandwiching(const Matrix &u, const Matrix &v, const Matrix &w, const Matrix &u1, const Matrix &v1, const Matrix &w1);
    void scale(int index, const Fraction &alpha, const Fraction &beta, const Fraction &gamma);
    void fixFractions();

    void copy(const FractionalScheme &scheme);
    void canonize();
    void saveJson(const std::string &path, bool withInvariants = false) const;
    void saveTxt(const std::string &path) const;
    void save(const std::string &path) const;
private:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const std::vector<Fraction> &u, const std::vector<Fraction> &v, const std::vector<Fraction> &w);

    bool validateEquation(int i, int j, int k) const;
    bool reconstructValue(int64_t a, int64_t mod, int64_t bound, Fraction &fraction) const;

    bool isEqualMatrices(int p, int index1, int index2) const;
    bool isInverseMatrices(int p, int index1, int index2) const;
    int compareMatrices(int p, int index1, int index2) const;
    bool isZeroMatrix(int p, int index) const;

    void selectFlip(FlipSet *flips, size_t index, int &i, int &j, int &k, int &index1, int &index2, std::mt19937 &generator);
    void flip(int i, int j, int k, int index1, int index2, bool inverse);
    void plus(int i, int j, int k, int index1, int index2, int variant);
    bool reduce(int i, int index1, int index2, int sign);
    void reduceAdd(int i, int index1, int index2);
    void reduceSub(int i, int index1, int index2);
    bool checkFlipReduce(int i, int j, int index1, int index2, int sign);

    std::vector<Fraction> addVectors(const std::vector<Fraction> &a, const std::vector<Fraction> &b) const;
    std::vector<Fraction> subVectors(const std::vector<Fraction> &a, const std::vector<Fraction> &b) const;

    int64_t gcdNumerators(const std::vector<Fraction> &fractions) const;
    int64_t lcmDenominators(const std::vector<Fraction> &fractions) const;

    void saveMatrix(std::ofstream &f, std::string name, const std::vector<Fraction> &values, int rows, int columns) const;
};
