#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_set>

#include "../algebra/fraction.h"
#include "base_scheme.h"

class FractionalScheme : public BaseScheme {
    std::vector<Fraction> uvw[3];
public:
    bool reconstruct(int n1, int n2, int n3, int rank, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, int64_t mod, int64_t bound);
    bool validate() const;

    bool isInteger() const;
    bool isTernary() const;

    int getComplexity() const;
    std::string getRing() const;
    std::string getUniqueValues() const;

    void canonize();
    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void save(const std::string &path) const;
private:
    bool validateEquation(int i, int j, int k) const;
    bool reconstructValue(int64_t a, int64_t mod, int64_t bound, Fraction &fraction) const;

    int64_t gcdNumerators(const std::vector<Fraction> &fractions) const;
    int64_t lcmDenominators(const std::vector<Fraction> &fractions) const;

    void saveMatrix(std::ofstream &f, std::string name, const std::vector<Fraction> &values, int rows, int columns) const;
};

bool FractionalScheme::reconstruct(int n1, int n2, int n3, int rank, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, int64_t mod, int64_t bound) {
    this->dimension[0] = n1;
    this->dimension[1] = n2;
    this->dimension[2] = n3;
    this->rank = rank;

    for (int i = 0; i < 3; i++) {
        elements[i] = dimension[i] * dimension[(i + 1) % 3];
        uvw[i].assign(rank * elements[i], 0);
    }

    for (int i = 0; i < rank * elements[0]; i++)
        if (!uvw[0][i].reconstruct(u[i], mod, bound))
            return false;

    for (int i = 0; i < rank * elements[1]; i++)
        if (!uvw[1][i].reconstruct(v[i], mod, bound))
            return false;

    for (int i = 0; i < rank * elements[2]; i++)
        if (!uvw[2][i].reconstruct(w[i], mod, bound))
            return false;

    return true;
}

bool FractionalScheme::validate() const {
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                if (!validateEquation(i, j, k))
                    return false;

    return true;
}

bool FractionalScheme::isInteger() const {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            if (!uvw[i][j].isInteger())
                return false;

    return true;
}

bool FractionalScheme::isTernary() const {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            if (!uvw[i][j].isTernaryInteger())
                return false;

    return true;
}

int FractionalScheme::getComplexity() const {
    int complexity = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            if (uvw[i][j] != 0)
                complexity++;

    return complexity - 2 * rank - elements[2];
}

std::string FractionalScheme::getRing() const {
    if (isTernary())
        return "ZT";

    if (isInteger())
        return "Z";

    return "Q";
}

std::string FractionalScheme::getUniqueValues() const {
    std::unordered_set<std::string> unique;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            unique.insert(uvw[i][j].pretty());

    std::vector<std::string> uniqueValues(unique.begin(), unique.end());
    std::sort(uniqueValues.begin(), uniqueValues.end());
    std::stringstream ss;

    ss << "{";
    for (size_t i = 0; i < uniqueValues.size(); i++)
        ss << (i > 0 ? ", " : "") << uniqueValues[i];
    ss << "}";

    return ss.str();
}

void FractionalScheme::canonize() {
    int64_t gcdV = gcdNumerators(uvw[1]);
    int64_t lcmV = lcmDenominators(uvw[1]);

    int64_t gcdW = gcdNumerators(uvw[2]);
    int64_t lcmW = lcmDenominators(uvw[2]);

    Fraction scaleV(gcdV, lcmV);
    Fraction scaleW(gcdW, lcmW);
    Fraction scaleU = scaleV * scaleW;

    for (auto &u : uvw[0])
        u *= scaleU;

    for (auto &v : uvw[1])
        v /= scaleV;

    for (auto &w : uvw[2])
        w /= scaleW;
}

void FractionalScheme::saveJson(const std::string &path) const {
    std::ofstream f(path);

    f << "{" << std::endl;
    f << "    \"n\": [" << dimension[0] << ", " << dimension[1] << ", " << dimension[2] << "]," << std::endl;
    f << "    \"m\": " << rank << "," << std::endl;
    f << "    \"z2\": false," << std::endl;
    f << "    \"complexity\": " << getComplexity() << "," << std::endl;

    saveMatrix(f, "u", uvw[0], rank, elements[0]);
    f << "," << std::endl;
    saveMatrix(f, "v", uvw[1], rank, elements[1]);
    f << "," << std::endl;
    saveMatrix(f, "w", uvw[2], rank, elements[2]);
    f << std::endl;
    f << "}" << std::endl;

    f.close();
}

void FractionalScheme::saveTxt(const std::string &path) const {
    std::ofstream f(path);
    f << dimension[0] << " " << dimension[1] << " " << dimension[2] << " " << rank << std::endl;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < rank * elements[i]; j++)
            f << (j > 0 ? " " : "") << uvw[i][j].numerator() << " " << uvw[i][j].denominator();

        f << std::endl;
    }

    f.close();
}

bool FractionalScheme::validateEquation(int i, int j, int k) const {
    int i1 = i / dimension[1];
    int i2 = i % dimension[1];
    int j1 = j / dimension[2];
    int j2 = j % dimension[2];
    int k1 = k / dimension[0];
    int k2 = k % dimension[0];

    Fraction target = (i2 == j1) && (i1 == k2) && (j2 == k1);
    Fraction equation = 0;

    for (int index = 0; index < rank; index++)
        equation += uvw[0][index * elements[0] + i] * uvw[1][index * elements[1] + j] * uvw[2][index * elements[2] + k];

    return equation == target;
}

int64_t FractionalScheme::gcdNumerators(const std::vector<Fraction> &fractions) const {
    int64_t result = 0;

    for (const auto &fraction: fractions) {
        int64_t num = fraction.numerator();

        if (num != 0)
            result = std::gcd(result, std::abs(num));
    }

    return result ? result : 1;
}

int64_t FractionalScheme::lcmDenominators(const std::vector<Fraction> &fractions) const {
    int64_t result = 1;

    for (const auto& fraction: fractions)
        result = std::lcm(result, fraction.denominator());

    return result;
}

void FractionalScheme::saveMatrix(std::ofstream &f, std::string name, const std::vector<Fraction> &values, int rows, int columns) const {
    f << "    \"" << name << "\": [" << std::endl;

    for (int i = 0; i < rows; i++) {
        f << "        [";

        for (int j = 0; j < columns; j++) {
            if (j > 0)
                f << ", ";

            if (values[i * columns + j].isInteger())
                f << values[i * columns + j].numerator();
            else
                f << '"' << values[i * columns + j] << '"';
        }

        f << "]" << (i < rows - 1 ? "," : "") << std::endl;
    }

    f << "    ]";
}
