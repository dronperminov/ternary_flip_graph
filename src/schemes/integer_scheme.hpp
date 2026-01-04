#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cassert>
#include <algorithm>

#include "../entities/integer_vector.hpp"
#include "../entities/flip_set.h"

template <typename T>
class IntegerScheme {
    int dimension[3];
    int elements[3];
    int rank;
    std::vector<IntegerVector<T>> uvw[3];
    FlipSet flips[3];
    std::vector<int> indices;

    std::uniform_int_distribution<int> boolDistribution;
    std::uniform_int_distribution<int> ijkDistribution;
public:
    IntegerScheme();
    IntegerScheme(const IntegerScheme<T> &scheme);

    bool initializeNaive(int n1, int n2, int n3);
    bool read(const std::string &path);

    int getRank() const;
    int getComplexity() const;
    int getDimension(int index) const;
    std::string getRing() const;
    int getAvailableFlips() const;

    bool tryFlip(std::mt19937 &generator);
    bool tryPlus(std::mt19937 &generator);
    bool trySplit(std::mt19937 &generator);
    bool tryExpand(std::mt19937 &generator);
    bool trySandwiching(std::mt19937 &generator) { return false; }
    bool tryReduce();

    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void copy(const IntegerScheme<T> &scheme);

    bool validate() const;
private:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const IntegerVector<T> &u, const IntegerVector<T> &v, const IntegerVector<T> &w);

    void flip(int i, int j, int k, int index1, int index2);
    bool plus(int i, int j, int k, int index1, int index2, int variant);
    void split(int i, int j, int k, int index1, int index2);
    void reduceAdd(int i, int index1, int index2);
    void reduceSub(int i, int index1, int index2);
    bool checkFlipReduce(int j, int k, int index1, int index2);

    bool fixSigns();
    bool validateDimensions() const;
    bool validateEquation(int i, int j, int k) const;
    void saveMatrix(std::ofstream &f, std::string name, const std::vector<IntegerVector<T>> &vectors) const;
};


template <typename T>
IntegerScheme<T>::IntegerScheme() : boolDistribution(0, 1), ijkDistribution(0, 2) {

}

template <typename T>
IntegerScheme<T>::IntegerScheme(const IntegerScheme<T> &scheme) {
    copy(scheme);
}

template <typename T>
bool IntegerScheme<T>::initializeNaive(int n1, int n2, int n3) {
    dimension[0] = n1;
    dimension[1] = n2;
    dimension[2] = n3;

    elements[0] = n1 * n2;
    elements[1] = n2 * n3;
    elements[2] = n3 * n1;

    rank = n1 * n2 * n3;

    if (!validateDimensions())
        return false;

    uvw[0].clear();
    uvw[1].clear();
    uvw[2].clear();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int k = 0; k < n2; k++) {
                uvw[0].emplace_back(IntegerVector<T>(n1 * n2, i * n2 + k));
                uvw[1].emplace_back(IntegerVector<T>(n2 * n3, k * n3 + j));
                uvw[2].emplace_back(IntegerVector<T>(n3 * n1, j * n1 + i));
            }
        }
    }

    initFlips();
    return true;
}

template <typename T>
bool IntegerScheme<T>::read(const std::string &path) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    f >> dimension[0] >> dimension[1] >> dimension[2] >> rank;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    if (!validateDimensions())
        return false;

    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++) {
            IntegerVector<T> vector(elements[i]);
            f >> vector;
            uvw[i].emplace_back(vector);
        }
    }

    f.close();

    if (!validate()) {
        std::cout << "Invalid scheme in the file \"" << path << "\"" << std::endl;
        return false;
    }

    initFlips();
    return true;
}

template <typename T>
int IntegerScheme<T>::getRank() const {
    return rank;
}

template <typename T>
int IntegerScheme<T>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += uvw[i][index].nonZeroCount();

    return count - 2 * rank - elements[2];
}

template <typename T>
int IntegerScheme<T>::getDimension(int index) const {
    return dimension[index];
}

template <typename T>
std::string IntegerScheme<T>::getRing() const {
    return "Z";
}

template <typename T>
int IntegerScheme<T>::getAvailableFlips() const {
    return flips[0].size() + flips[1].size() + flips[2].size();
}

template <typename T>
bool IntegerScheme<T>::tryFlip(std::mt19937 &generator) {
    size_t size = flips[0].size() + flips[1].size() + flips[2].size();

    if (!size)
        return false;

    indices.resize(size);
    for (size_t i = 0; i < size; i++)
        indices[i] = i;

    for (size_t p = 0; p < size; p++) {
        size_t q = generator() % (size - p);
        size_t index = indices[q];
        std::swap(indices[p], indices[q]);

        int i, j, k;

        if (index < flips[0].size()) {
            i = 0;
            j = 1;
            k = 2;
        }
        else if (index < flips[0].size() + flips[1].size()) {
            i = 1;
            j = 0;
            k = 2;
            index -= flips[0].size();
        }
        else {
            i = 2;
            j = 0;
            k = 1;
            index -= flips[0].size() + flips[1].size();
        }

        int index1 = flips[i].index1(index);
        int index2 = flips[i].index2(index);

        if (boolDistribution(generator))
            std::swap(j, k);

        if (boolDistribution(generator))
            std::swap(index1, index2);

        if (uvw[j][index1].limitSum(uvw[j][index2], j != 2) && uvw[k][index2].limitSub(uvw[k][index1], false)) {
            if (k == 2 || uvw[k][index2].positiveFirstNonZeroSub(uvw[k][index1]))
                flip(i, j, k, index1, index2);
            else
                flip(i, j, k, index2, index1);
            return true;
        }

        if (uvw[k][index1].limitSum(uvw[k][index2], k != 2) && uvw[j][index2].limitSub(uvw[j][index1], false)) {
            if (j == 2 || uvw[j][index2].positiveFirstNonZeroSub(uvw[j][index1]))
                flip(i, k, j, index1, index2);
            else
                flip(i, k, j, index2, index1);
            return true;
        }
    }

    return false;
}

template <typename T>
bool IntegerScheme<T>::tryPlus(std::mt19937 &generator) {
    std::uniform_int_distribution<int> distribution(0, rank - 1);
    int index1 = distribution(generator);
    int index2 = distribution(generator);

    while (index1 == index2 || uvw[0][index1] == uvw[0][index2] || uvw[1][index1] == uvw[1][index2] || uvw[2][index1] == uvw[2][index2]) {
        index1 = distribution(generator);
        index2 = distribution(generator);
    }

    int permutation[3] = {0, 1, 2};
    std::shuffle(permutation, permutation + 3, generator);

    return plus(permutation[0], permutation[1], permutation[2], index1, index2, ijkDistribution(generator));
}

template <typename T>
bool IntegerScheme<T>::trySplit(std::mt19937 &generator) {
    std::uniform_int_distribution<int> distribution(0, rank - 1);
    int index1, index2, i;

    do {
        index1 = distribution(generator);
        index2 = distribution(generator);
        i = ijkDistribution(generator);
    } while (index1 == index2 || uvw[i][index1] == uvw[i][index2]);

    if (!uvw[i][index1].limitSub(uvw[i][index2], false))
        return false;

    if (i == 2 || uvw[i][index1].positiveFirstNonZeroSub(uvw[i][index2]))
        split(i, (i + 1) % 3, (i + 2) % 3, index1, index2);
    else
        split(i, (i + 1) % 3, (i + 2) % 3, index2, index1);

    return true;
}

template <typename T>
bool IntegerScheme<T>::tryExpand(std::mt19937 &generator) {
    if (rank >= dimension[0] * dimension[1] * dimension[2])
        return false;

    if (boolDistribution(generator))
        return tryPlus(generator);

    return trySplit(generator);
}

template <typename T>
bool IntegerScheme<T>::tryReduce() {
    for (size_t i = 0; i < flips[0].size(); i++) {
        int index1 = flips[0].index1(i);
        int index2 = flips[0].index2(i);

        if (uvw[1][index1] == uvw[1][index2] && uvw[2][index1].limitSum(uvw[2][index2], false)) {
            reduceAdd(2, index1, index2);
            return true;
        }

        int cmp2 = uvw[2][index1].compare(uvw[2][index2]);
        if (cmp2 == 1 && uvw[1][index1].limitSum(uvw[1][index2], true)) {
            reduceAdd(1, index1, index2);
            return true;
        }

        if (cmp2 == -1 && uvw[1][index1].limitSub(uvw[1][index2], true)) {
            reduceSub(1, index1, index2);
            return true;
        }
    }

    for (size_t i = 0; i < flips[1].size(); i++) {
        int index1 = flips[1].index1(i);
        int index2 = flips[1].index2(i);
        int cmp2 = uvw[2][index1].compare(uvw[2][index2]);

        if (cmp2 == 1 && uvw[0][index1].limitSum(uvw[0][index2], true)) {
            reduceAdd(0, index1, index2);
            return true;
        }

        if (cmp2 == -1 && uvw[0][index1].limitSub(uvw[0][index2], true)) {
            reduceSub(0, index1, index2);
            return true;
        }
    }

    return false;
}

template <typename T>
bool IntegerScheme<T>::validate() const {
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                if (!validateEquation(i, j, k))
                    return false;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            if (!uvw[i][index].limit(false))
                return false;

    return true;
}

template <typename T>
void IntegerScheme<T>::saveJson(const std::string &path) const {
    std::ofstream f(path);

    f << "{" << std::endl;
    f << "    \"n\": [" << dimension[0] << ", " << dimension[1] << ", " << dimension[2] << "]," << std::endl;
    f << "    \"m\": " << rank << "," << std::endl;
    f << "    \"z2\": false," << std::endl;
    f << "    \"complexity\": " << getComplexity() << "," << std::endl;

    saveMatrix(f, "u", uvw[0]);
    f << "," << std::endl;
    saveMatrix(f, "v", uvw[1]);
    f << "," << std::endl;
    saveMatrix(f, "w", uvw[2]);
    f << std::endl;
    f << "}" << std::endl;

    f.close();
}

template <typename T>
void IntegerScheme<T>::saveTxt(const std::string &path) const {
    std::ofstream f(path);

    f << dimension[0] << " " << dimension[1] << " " << dimension[2] << " " << rank << std::endl;
    
    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++)
            for (int j = 0; j < elements[i]; j++)
                f << uvw[i][index][j] << " ";

        f << std::endl;
    }

    f.close();
}

template <typename T>
void IntegerScheme<T>::copy(const IntegerScheme<T> &scheme) {
    rank = scheme.rank;

    for (int i = 0; i < 3; i++) {
        dimension[i] = scheme.dimension[i];
        elements[i] = scheme.elements[i];
        uvw[i].clear();

        for (int index = 0; index < rank; index++) {
            IntegerVector<T> vector(elements[i]);

            for (int j = 0; j < elements[i]; j++)
                vector.set(j, scheme.uvw[i][index][j]);

            uvw[i].emplace_back(vector);
        }
    }

    initFlips();
}

template <typename T>
void IntegerScheme<T>::initFlips() {
    for (int i = 0; i < 3; i++) {
        flips[i].clear();

        for (int index1 = 0; index1 < rank; index1++)
            for (int index2 = index1 + 1; index2 < rank; index2++)
                if (uvw[i][index1] == uvw[i][index2])
                    flips[i].add(index1, index2);
    }
}

template <typename T>
void IntegerScheme<T>::removeZeroes() {
    for (int index = 0; index < rank; index++)
        if (!uvw[0][index] || !uvw[1][index] || !uvw[2][index])
            removeAt(index--);
}

template <typename T>
void IntegerScheme<T>::removeAt(int index) {
    if (index != rank) {
        uvw[0][index] = uvw[0].back();
        uvw[1][index] = uvw[1].back();
        uvw[2][index] = uvw[2].back();
    }

    rank--;
    uvw[0].pop_back();
    uvw[1].pop_back();
    uvw[2].pop_back();
}

template <typename T>
void IntegerScheme<T>::addTriplet(int i, int j, int k, const IntegerVector<T> &u, const IntegerVector<T> &v, const IntegerVector<T> &w) {
    uvw[i].emplace_back(u);
    uvw[j].emplace_back(v);
    uvw[k].emplace_back(w);
    rank++;
}

template <typename T>
void IntegerScheme<T>::flip(int i, int j, int k, int index1, int index2) {
    uvw[j][index1] += uvw[j][index2];
    uvw[k][index2] -= uvw[k][index1];

    flips[j].remove(index1);
    flips[k].remove(index2);

    if (!uvw[j][index1] || !uvw[k][index2]) {
        removeZeroes();
        initFlips();

        while (tryReduce())
            ;
        return;
    }

    for (int index = 0; index < rank; index++) {
        if (index != index1 && uvw[j][index] == uvw[j][index1]) {
            if (checkFlipReduce(i, k, index, index1))
                return;

            flips[j].add(index1, index);
        }

        if (index != index2 && uvw[k][index] == uvw[k][index2]) {
            if (checkFlipReduce(i, j, index, index2))
                return;

            flips[k].add(index2, index);
        }
    }
}

template <typename T>
bool IntegerScheme<T>::plus(int i, int j, int k, int index1, int index2, int variant) {
    IntegerVector<T> a1(uvw[i][index1]);
    IntegerVector<T> b1(uvw[j][index1]);
    IntegerVector<T> c1(uvw[k][index1]);

    IntegerVector<T> a2(uvw[i][index2]);
    IntegerVector<T> b2(uvw[j][index2]);
    IntegerVector<T> c2(uvw[k][index2]);

    IntegerVector<T> aAdd = a1 + a2;
    IntegerVector<T> bAdd = b1 + b2;
    IntegerVector<T> cAdd = c1 + c2;

    IntegerVector<T> aSub = a2 - a1;
    IntegerVector<T> bSub = b2 - b1;
    IntegerVector<T> cSub = c2 - c1;

    if (variant == 0 && aSub.limit(i != 2) && bAdd.limit(j != 2) && cSub.limit(k != 2)) {
        uvw[j][index1] = bAdd;
        uvw[i][index2] = aSub;
        addTriplet(i, j, k, a1, b2, cSub);
    }
    else if (variant == 1 && aSub.limit(i != 2) && bSub.limit(j != 2) && cAdd.limit(k != 2)) {
        uvw[k][index1] = cAdd;
        uvw[j][index2] = bSub;
        addTriplet(i, j, k, aSub, b1, c2);
    }
    else if (aAdd.limit(i != 2) && bSub.limit(j != 2) && cSub.limit(k != 2)) {
        uvw[i][index1] = aAdd;
        uvw[k][index2] = cSub;
        addTriplet(i, j, k, a2, bSub, c1);
    }
    else
        return false;

    removeZeroes();
    fixSigns();
    initFlips();
    return true;
}

template <typename T>
void IntegerScheme<T>::split(int i, int j, int k, int index1, int index2) {
    IntegerVector<T> u = uvw[i][index1] - uvw[i][index2];
    IntegerVector<T> v(uvw[j][index1]);
    IntegerVector<T> w(uvw[k][index1]);

    addTriplet(i, j, k, u, v, w);
    uvw[i][index1] = uvw[i][index2];

    removeZeroes();
    fixSigns();
    initFlips();
}

template <typename T>
void IntegerScheme<T>::reduceAdd(int i, int index1, int index2) {
    uvw[i][index1] += uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <typename T>
void IntegerScheme<T>::reduceSub(int i, int index1, int index2) {
    uvw[i][index1] -= uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <typename T>
bool IntegerScheme<T>::checkFlipReduce(int i, int j, int index1, int index2) {
    int cmpI = uvw[i][index1].compare(uvw[i][index2]);
    if (cmpI == 1 && uvw[j][index1].limitSum(uvw[j][index2], j != 2)) {
        reduceAdd(j, index1, index2);
        return true;
    }

    if (cmpI == -1 && uvw[j][index1].limitSub(uvw[j][index2], false)) {
        if (j == 2 || uvw[j][index1].positiveFirstNonZeroSub(uvw[j][index2]))
            reduceSub(j, index1, index2);
        else
            reduceSub(j, index2, index1);
        return true;
    }

    int cmpJ = uvw[j][index1].compare(uvw[j][index2]);
    if (cmpJ == 1 && uvw[i][index1].limitSum(uvw[i][index2], i != 2)) {
        reduceAdd(i, index1, index2);
        return true;
    }

    if (cmpJ == -1 && uvw[i][index1].limitSub(uvw[i][index2], false)) {
        if (i == 2 || uvw[i][index1].positiveFirstNonZeroSub(uvw[i][index2]))
            reduceSub(i, index1, index2);
        else
            reduceSub(i, index2, index1);
        return true;
    }

    return false;
}

template <typename T>
bool IntegerScheme<T>::fixSigns() {
    bool changed = false;

    for (int index = 0; index < rank; index++) {
        bool i = uvw[0][index].positiveFirstNonZero();
        bool j = uvw[1][index].positiveFirstNonZero();

        if (i && j)
            continue;

        if (!i && !j) {
            uvw[0][index].inverse();
            uvw[1][index].inverse();
        }
        else if (!i) {
            uvw[0][index].inverse();
            uvw[2][index].inverse();
        }
        else {
            uvw[1][index].inverse();
            uvw[2][index].inverse();
        }

        changed = true;
    }

    return changed;
}

template <typename T>
bool IntegerScheme<T>::validateDimensions() const {
    for (int i = 0; i < 3; i++) {
        if (dimension[i] < 1) {
            std::cout << "Invalid dimension \"" << dimension[i] << "\"." << std::endl;
            return false;
        }
    }

    if (rank < 1) {
        std::cout << "Invalid rank \"" << rank << "\"" << std::endl;
        return false;
    }

    return true;
}

template <typename T>
bool IntegerScheme<T>::validateEquation(int i, int j, int k) const {
    int i1 = i / dimension[1];
    int i2 = i % dimension[1];
    int j1 = j / dimension[2];
    int j2 = j % dimension[2];
    int k1 = k / dimension[0];
    int k2 = k % dimension[0];

    int target = (i2 == j1) && (i1 == k2) && (j2 == k1);
    int equation = 0;

    for (int index = 0; index < rank; index++)
        equation += uvw[0][index][i] * uvw[1][index][j] * uvw[2][index][k];

    return equation == target;
}

template <typename T>
void IntegerScheme<T>::saveMatrix(std::ofstream &f, std::string name, const std::vector<IntegerVector<T>> &vectors) const {
    f << "    \"" << name << "\": [" << std::endl;

    for (size_t index = 0; index < vectors.size(); index++)
        f << "        [" << vectors[index] << "]" << (index < vectors.size() - 1 ? "," : "") << std::endl;

    f << "    ]";
}
