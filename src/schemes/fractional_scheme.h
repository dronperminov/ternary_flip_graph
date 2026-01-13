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
#include "base_scheme.h"

class FractionalScheme : public BaseScheme {
    std::vector<Fraction> uvw[3];
public:
    bool reconstruct(int n1, int n2, int n3, int rank, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, int64_t mod, int64_t bound);
    bool validate() const;

    bool isInteger() const;
    bool isTernary() const;

    int getFractionsCount() const;
    int getComplexity() const;
    int64_t getWeight() const;
    std::string getRing() const;
    std::string getUniqueValues() const;

    std::string getTypeInvariant() const;

    void copy(const FractionalScheme &scheme);
    void canonize();
    void saveJson(const std::string &path, bool withInvariants = false) const;
    void saveTxt(const std::string &path) const;
    void save(const std::string &path) const;
private:
    bool validateEquation(int i, int j, int k) const;
    bool reconstructValue(int64_t a, int64_t mod, int64_t bound, Fraction &fraction) const;

    int64_t gcdNumerators(const std::vector<Fraction> &fractions) const;
    int64_t lcmDenominators(const std::vector<Fraction> &fractions) const;

    void saveMatrix(std::ofstream &f, std::string name, const std::vector<Fraction> &values, int rows, int columns) const;
};
