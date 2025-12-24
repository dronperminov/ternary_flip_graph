#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cassert>
#include <algorithm>

#include "../entities/ternary_vector.hpp"
#include "base_scheme.h"

template <typename T>
class TernaryScheme : public BaseScheme {
    std::vector<TernaryVector<T>> uvw[3];
    std::vector<int> indices;
public:
    TernaryScheme();
    TernaryScheme(const TernaryScheme<T> &scheme);

    bool initializeNaive(int n1, int n2, int n3);
    bool read(const std::string &path);
    bool read(std::istream &is);

    int getComplexity() const;
    std::string getRing() const;
    std::string getHash() const;

    bool tryFlip(std::mt19937 &generator);
    bool tryPlus(std::mt19937 &generator);
    bool trySplit(std::mt19937 &generator);
    bool tryExpand(std::mt19937 &generator);
    bool tryReduce();

    bool tryProject(std::mt19937 &generator, int minN);
    bool tryExtend(std::mt19937 &generator, int maxN, int maxRank);
    bool tryMerge(const TernaryScheme<T> &scheme, std::mt19937 &generator, int maxN, int maxRank);
    bool tryProduct(const TernaryScheme<T> &scheme, int maxN, int maxRank);

    void swapSizes(std::mt19937 &generator);
    void swapSizes(int p1, int p2);
    void merge(const TernaryScheme<T> &scheme, int p);
    void project(int p, int q);
    void extend(int p);
    void product(const TernaryScheme<T> &scheme);

    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void copy(const TernaryScheme<T> &scheme);

    bool validate() const;
private:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const TernaryVector<T> &u, const TernaryVector<T> &v, const TernaryVector<T> &w);
    void excludeColumn(int matrix, int column);
    void excludeRow(int matrix, int row);
    void addColumn(int matrix);
    void addRow(int matrix);

    void flip(int i, int j, int k, int index1, int index2);
    bool plus(int i, int j, int k, int index1, int index2, int variant);
    void split(int i, int j, int k, int index1, int index2);
    void reduceAdd(int i, int index1, int index2);
    void reduceSub(int i, int index1, int index2);
    bool checkFlipReduce(int j, int k, int index1, int index2);

    bool isValidProject(int p, int minN) const;
    bool isValidExtension(int p, int maxN, int maxRank) const;
    bool isValidProduct(const TernaryScheme<T> &scheme, int maxN, int maxRank) const;
    bool isValidMerge(int p, const TernaryScheme<T> &scheme, int maxN, int maxRank) const;

    bool fixSigns();
    bool validateDimensions() const;
    bool validateEquation(int i, int j, int k) const;
    void saveMatrix(std::ofstream &f, std::string name, const std::vector<TernaryVector<T>> &vectors) const;
};

template <typename T>
TernaryScheme<T>::TernaryScheme() {

}

template <typename T>
TernaryScheme<T>::TernaryScheme(const TernaryScheme<T> &scheme) {
    copy(scheme);
}

template <typename T>
bool TernaryScheme<T>::initializeNaive(int n1, int n2, int n3) {
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
                uvw[0].emplace_back(TernaryVector<T>(n1 * n2, i * n2 + k));
                uvw[1].emplace_back(TernaryVector<T>(n2 * n3, k * n3 + j));
                uvw[2].emplace_back(TernaryVector<T>(n3 * n1, j * n1 + i));
            }
        }
    }

    initFlips();
    return true;
}

template <typename T>
bool TernaryScheme<T>::read(const std::string &path) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    bool valid = read(f);
    f.close();

    if (!valid) {
        std::cout << "Invalid scheme in the file \"" << path << "\"" << std::endl;
        return false;
    }

    return true;
}

template <typename T>
bool TernaryScheme<T>::read(std::istream &is) {
    is >> dimension[0] >> dimension[1] >> dimension[2] >> rank;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    if (!validateDimensions())
        return false;

    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++) {
            TernaryVector<T> vector(elements[i]);
            is >> vector;
            uvw[i].emplace_back(vector);
        }
    }

    if (!validate())
        return false;

    fixSigns();
    initFlips();
    return true;
}

template <typename T>
int TernaryScheme<T>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += uvw[i][index].nonZeroCount();

    return count - 2 * rank - elements[2];
}

template <typename T>
std::string TernaryScheme<T>::getRing() const {
    return "ZT";
}

template <typename T>
std::string TernaryScheme<T>::getHash() const {
    std::vector<std::string> lines;

    for (int index = 0; index < rank; index++) {
        std::stringstream ss;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < elements[i]; j++)
                ss << (uvw[i][index][j] + 1);

        lines.push_back(ss.str());
    }

    std::sort(lines.begin(), lines.end());
    std::stringstream hash;
    for (int index = 0; index < rank; index++)
        hash << lines[index];

    return hash.str();
}

template <typename T>
bool TernaryScheme<T>::tryFlip(std::mt19937 &generator) {
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
bool TernaryScheme<T>::tryPlus(std::mt19937 &generator) {
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
bool TernaryScheme<T>::trySplit(std::mt19937 &generator) {
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
bool TernaryScheme<T>::tryExpand(std::mt19937 &generator) {
    if (rank >= dimension[0] * dimension[1] * dimension[2])
        return false;

    if (boolDistribution(generator))
        return tryPlus(generator);

    return trySplit(generator);
}

template <typename T>
bool TernaryScheme<T>::tryReduce() {
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
bool TernaryScheme<T>::tryProject(std::mt19937 &generator, int minN) {
    std::vector<int> indices;

    for (int i = 0; i < 3; i++)
        if (isValidProject(i, minN))
            indices.push_back(i);

    if (!indices.size())
        return false;

    int p = indices[generator() % indices.size()];
    int q = generator() % dimension[p];
    project(p, q);

    while (tryReduce())
        ;

    return true;
}

template <typename T>
bool TernaryScheme<T>::tryExtend(std::mt19937 &generator, int maxN, int maxRank) {
    std::vector<int> indices;

    for (int i = 0; i < 3; i++)
        if (isValidExtension(i, maxN, maxRank))
            indices.push_back(i);

    if (!indices.size())
        return false;

    extend(indices[generator() % indices.size()]);
    return true;
}

template <typename T>
bool TernaryScheme<T>::tryMerge(const TernaryScheme<T> &scheme, std::mt19937 &generator, int maxN, int maxRank) {
    std::vector<int> indices;

    for (int i = 0; i < 3; i++)
        if (isValidMerge(i, scheme, maxN, maxRank))
            indices.push_back(i);

    if (!indices.size())
        return false;

    merge(scheme, indices[generator() % indices.size()]);
    return true;
}

template <typename T>
bool TernaryScheme<T>::tryProduct(const TernaryScheme<T> &scheme, int maxN, int maxRank) {
    if (!isValidProduct(scheme, maxN, maxRank))
        return false;

    product(scheme);
    return true;
}

template <typename T>
void TernaryScheme<T>::swapSizes(std::mt19937 &generator) {
    int p1, p2;

    do {
        p1 = ijkDistribution(generator);
        p2 = ijkDistribution(generator);
    } while (p1 == p2);

    swapSizes(p1, p2);
}

template <typename T>
void TernaryScheme<T>::swapSizes(int p1, int p2) {
    if (p1 == p2)
        return;

    if (p1 > p2)
        std::swap(p1, p2);

    int indices[3] = {2, 0, 1};
    int dimensionNew[3];

    std::swap(indices[p1], indices[p2]);

    for (int i = 0; i < 3; i++)
        dimensionNew[i] = dimension[(indices[i] + 1) % 3];

    for (int index = 0; index < rank; index++) {
        TernaryVector<T> u(dimensionNew[0] * dimensionNew[1]);
        TernaryVector<T> v(dimensionNew[1] * dimensionNew[2]);
        TernaryVector<T> w(dimensionNew[2] * dimensionNew[0]);

        for (int i = 0; i < dimensionNew[0]; i++)
            for (int j = 0; j < dimensionNew[1]; j++)
                u.set(i * dimensionNew[1] + j, uvw[indices[0]][index][j * dimensionNew[0] + i]);

        for (int i = 0; i < dimensionNew[1]; i++)
            for (int j = 0; j < dimensionNew[2]; j++)
                v.set(i * dimensionNew[2] + j, uvw[indices[1]][index][j * dimensionNew[1] + i]);

        for (int i = 0; i < dimensionNew[2]; i++)
            for (int j = 0; j < dimensionNew[0]; j++)
                w.set(i * dimensionNew[0] + j, uvw[indices[2]][index][j * dimensionNew[2] + i]);

        uvw[0][index] = u;
        uvw[1][index] = v;
        uvw[2][index] = w;
    }

    for (int i = 0; i < 3; i++) {
        dimension[i] = dimensionNew[i];
        elements[i] = dimensionNew[i] * dimensionNew[(i + 1) % 3];
    }

    fixSigns();
    initFlips();
}

template <typename T>
void TernaryScheme<T>::merge(const TernaryScheme<T> &scheme, int p) {
    int dimensionNew[3];
    int elementsNew[3];
    int d[3];

    for (int i = 0; i < 3; i++) {
        dimensionNew[i] = i == p ? dimension[i] + scheme.dimension[i] : dimension[i];
        d[i] = i == p ? dimension[i] : 0;
    }

    for (int i = 0; i < 3; i++)
        elementsNew[i] = dimensionNew[i] * dimensionNew[(i + 1) % 3];

    for (int index = 0; index < rank; index++) {
        TernaryVector<T> u(elementsNew[0]);
        TernaryVector<T> v(elementsNew[1]);
        TernaryVector<T> w(elementsNew[2]);

        for (int i = 0; i < dimension[0]; i++)
            for (int j = 0; j < dimension[1]; j++)
                u.set(i * dimensionNew[1] + j, uvw[0][index][i * dimension[1] + j]);

        for (int i = 0; i < dimension[1]; i++)
            for (int j = 0; j < dimension[2]; j++)
                v.set(i * dimensionNew[2] + j, uvw[1][index][i * dimension[2] + j]);

        for (int i = 0; i < dimension[2]; i++)
            for (int j = 0; j < dimension[0]; j++)
                w.set(i * dimensionNew[0] + j, uvw[2][index][i * dimension[0] + j]);

        uvw[0][index] = u;
        uvw[1][index] = v;
        uvw[2][index] = w;
    }

    for (int index = 0; index < scheme.rank; index++) {
        TernaryVector<T> u(elementsNew[0]);
        TernaryVector<T> v(elementsNew[1]);
        TernaryVector<T> w(elementsNew[2]);

        for (int i = 0; i < scheme.dimension[0]; i++)
            for (int j = 0; j < scheme.dimension[1]; j++)
                u.set((i + d[0]) * dimensionNew[1] + j + d[1], scheme.uvw[0][index][i * scheme.dimension[1] + j]);

        for (int i = 0; i < scheme.dimension[1]; i++)
            for (int j = 0; j < scheme.dimension[2]; j++)
                v.set((i + d[1]) * dimensionNew[2] + j + d[2], scheme.uvw[1][index][i * scheme.dimension[2] + j]);

        for (int i = 0; i < scheme.dimension[2]; i++)
            for (int j = 0; j < scheme.dimension[0]; j++)
                w.set((i + d[2]) * dimensionNew[0] + j + d[0], scheme.uvw[2][index][i * scheme.dimension[0] + j]);

        addTriplet(0, 1, 2, u, v, w);
    }

    for (int i = 0; i < 3; i++) {
        dimension[i] = dimensionNew[i];
        elements[i] = elementsNew[i];
    }

    initFlips();
}

template <typename T>
void TernaryScheme<T>::project(int p, int q) {
    excludeRow(p, q);
    excludeColumn((p + 2) % 3, q);
    dimension[p]--;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    removeZeroes();
    fixSigns();
    initFlips();
}

template <typename T>
void TernaryScheme<T>::extend(int p) {
    addRow(p);
    addColumn((p + 2) % 3);

    if (p == 0) {
        for (int i = 0; i < dimension[2]; i++) {
            for (int j = 0; j < dimension[1]; j++) {
                TernaryVector<T> u((dimension[0] + 1) * dimension[1], dimension[0] * dimension[1] + j);
                TernaryVector<T> v(dimension[1] * dimension[2], j * dimension[2] + i);
                TernaryVector<T> w(dimension[2] * (dimension[0] + 1), i * (dimension[0] + 1) + dimension[0]);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }
    else if (p == 1) {
        for (int i = 0; i < dimension[0]; i++) {
            for (int j = 0; j < dimension[2]; j++) {
                TernaryVector<T> u(dimension[0] * (dimension[1] + 1), i * (dimension[1] + 1) + dimension[1]);
                TernaryVector<T> v((dimension[1] + 1) * dimension[2], dimension[1] * dimension[2] + j);
                TernaryVector<T> w(dimension[2] * dimension[0], j * dimension[0] + i);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }
    else {
        for (int i = 0; i < dimension[0]; i++) {
            for (int j = 0; j < dimension[1]; j++) {
                TernaryVector<T> u(dimension[0] * dimension[1], i * dimension[1] + j);
                TernaryVector<T> v(dimension[1] * (dimension[2] + 1), j * (dimension[2] + 1) + dimension[2]);
                TernaryVector<T> w((dimension[2] + 1) * dimension[0], dimension[2] * dimension[0] + i);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }

    dimension[p]++;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    initFlips();
}

template <typename T>
void TernaryScheme<T>::product(const TernaryScheme<T> &scheme2) {
    TernaryScheme<T> scheme1(*this);

    for (int i = 0; i < 3; i++)
        dimension[i] = scheme1.dimension[i] * scheme2.dimension[i];

    for (int i = 0; i < 3; i++) {
        elements[i] = dimension[i] * dimension[(i + 1) % 3];
        uvw[i].clear();
    }

    rank = scheme1.rank * scheme2.rank;

    for (int index1 = 0; index1 < scheme1.rank; index1++) {
        for (int index2 = 0; index2 < scheme2.rank; index2++) {
            for (int p = 0; p < 3; p++) {
                int p1 = (p + 1) % 3;

                TernaryVector<T> vector(elements[p]);

                for (int i = 0; i < scheme1.elements[p]; i++) {
                    for (int j = 0; j < scheme2.elements[p]; j++) {
                        int row1 = i / scheme1.dimension[p1];
                        int col1 = i % scheme1.dimension[p1];
                        int value1 = scheme1.uvw[p][index1][i];

                        int row2 = j / scheme2.dimension[p1];
                        int col2 = j % scheme2.dimension[p1];
                        int value2 = scheme2.uvw[p][index2][j];

                        int row = row1 * scheme2.dimension[p] + row2;
                        int col = col1 * scheme2.dimension[p1] + col2;

                        vector.set(row * dimension[p1] + col, value1 * value2);
                    }
                }

                uvw[p].push_back(vector);
            }
        }
    }

    initFlips();
}

template <typename T>
void TernaryScheme<T>::saveJson(const std::string &path) const {
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
void TernaryScheme<T>::saveTxt(const std::string &path) const {
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
void TernaryScheme<T>::copy(const TernaryScheme<T> &scheme) {
    rank = scheme.rank;

    for (int i = 0; i < 3; i++) {
        dimension[i] = scheme.dimension[i];
        elements[i] = scheme.elements[i];
        uvw[i].clear();

        for (int index = 0; index < rank; index++) {
            TernaryVector<T> vector(elements[i]);

            for (int j = 0; j < elements[i]; j++)
                vector.set(j, scheme.uvw[i][index][j]);

            uvw[i].emplace_back(vector);
        }
    }

    initFlips();
}

template <typename T>
bool TernaryScheme<T>::validate() const {
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
void TernaryScheme<T>::initFlips() {
    for (int i = 0; i < 3; i++) {
        flips[i].clear();

        for (int index1 = 0; index1 < rank; index1++)
            for (int index2 = index1 + 1; index2 < rank; index2++)
                if (uvw[i][index1] == uvw[i][index2])
                    flips[i].add(index1, index2);
    }
}

template <typename T>
void TernaryScheme<T>::removeZeroes() {
    for (int index = 0; index < rank; index++)
        if (!uvw[0][index] || !uvw[1][index] || !uvw[2][index])
            removeAt(index--);
}

template <typename T>
void TernaryScheme<T>::removeAt(int index) {
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
void TernaryScheme<T>::addTriplet(int i, int j, int k, const TernaryVector<T> &u, const TernaryVector<T> &v, const TernaryVector<T> &w) {
    uvw[i].emplace_back(u);
    uvw[j].emplace_back(v);
    uvw[k].emplace_back(w);
    rank++;
}

template <typename T>
void TernaryScheme<T>::excludeColumn(int matrix, int column) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];
    std::vector<int> oldColumns;

    for (int j = 0; j < n2; j++)
        if (j != column)
            oldColumns.push_back(j);

    for (int index = 0; index < rank; index++) {
        TernaryVector<T> vector(n1 * (n2 - 1));

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2 - 1; j++)
                vector.set(i * (n2 - 1) + j, uvw[matrix][index][i * n2 + oldColumns[j]]);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void TernaryScheme<T>::excludeRow(int matrix, int row) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];
    std::vector<int> oldRows;

    for (int i = 0; i < n1; i++)
        if (i != row)
            oldRows.push_back(i);

    for (int index = 0; index < rank; index++) {
        TernaryVector<T> vector((n1 - 1) * n2);

        for (int i = 0; i < n1 - 1; i++)
            for (int j = 0; j < n2; j++)
                vector.set(i * n2 + j, uvw[matrix][index][oldRows[i] * n2 + j]);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void TernaryScheme<T>::addColumn(int matrix) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];

    for (int index = 0; index < rank; index++) {
        TernaryVector<T> vector(n1 * (n2 + 1));

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vector.set(i * (n2 + 1) + j, uvw[matrix][index][i * n2 + j]);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void TernaryScheme<T>::addRow(int matrix) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];

    for (int index = 0; index < rank; index++) {
        TernaryVector<T> vector((n1 + 1) * n2);

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vector.set(i * n2 + j, uvw[matrix][index][i * n2 + j]);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void TernaryScheme<T>::flip(int i, int j, int k, int index1, int index2) {
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
bool TernaryScheme<T>::plus(int i, int j, int k, int index1, int index2, int variant) {
    TernaryVector<T> a1(uvw[i][index1]);
    TernaryVector<T> b1(uvw[j][index1]);
    TernaryVector<T> c1(uvw[k][index1]);

    TernaryVector<T> a2(uvw[i][index2]);
    TernaryVector<T> b2(uvw[j][index2]);
    TernaryVector<T> c2(uvw[k][index2]);

    TernaryVector<T> aAdd = a1 + a2;
    TernaryVector<T> bAdd = b1 + b2;
    TernaryVector<T> cAdd = c1 + c2;

    TernaryVector<T> aSub = a2 - a1;
    TernaryVector<T> bSub = b2 - b1;
    TernaryVector<T> cSub = c2 - c1;

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
void TernaryScheme<T>::split(int i, int j, int k, int index1, int index2) {
    TernaryVector<T> u = uvw[i][index1] - uvw[i][index2];
    TernaryVector<T> v(uvw[j][index1]);
    TernaryVector<T> w(uvw[k][index1]);

    addTriplet(i, j, k, u, v, w);
    uvw[i][index1] = uvw[i][index2];

    removeZeroes();
    fixSigns();
    initFlips();
}

template <typename T>
void TernaryScheme<T>::reduceAdd(int i, int index1, int index2) {
    uvw[i][index1] += uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <typename T>
void TernaryScheme<T>::reduceSub(int i, int index1, int index2) {
    uvw[i][index1] -= uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <typename T>
bool TernaryScheme<T>::checkFlipReduce(int i, int j, int index1, int index2) {
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
bool TernaryScheme<T>::isValidProject(int p, int minN) const {
    return dimension[p] > minN && dimension[(p + 1) % 3] >= minN && dimension[(p + 2) % 3] >= minN;
}

template <typename T>
bool TernaryScheme<T>::isValidExtension(int p, int maxN, int maxRank) const {
    if (rank + dimension[(p + 1) % 3] * dimension[(p + 2) % 3] > maxRank)
        return false;

    int dimensionNew[3] = {dimension[0], dimension[1], dimension[2]};
    int maxElements = sizeof(T) * 8;
    dimensionNew[p]++;

    for (int i = 0; i < 3; i++) {
        if (dimensionNew[i] * dimensionNew[(i + 1) % 3] > maxElements)
            return false;

        if (dimensionNew[i] > maxN)
            return false;
    }

    return true;
}

template <typename T>
bool TernaryScheme<T>::isValidProduct(const TernaryScheme<T> &scheme, int maxN, int maxRank) const {
    if (rank * scheme.rank > maxRank)
        return false;

    int dimensionNew[3];

    for (int i = 0; i < 3; i++) {
        dimensionNew[i] = dimension[i] * scheme.dimension[i];

        if (dimensionNew[i] > maxN)
            return false;
    }

    int maxElements = sizeof(T) * 8;
    for (int i = 0; i < 3; i++)
        if (dimensionNew[i] * dimensionNew[(i + 1) % 3] > maxElements)
            return false;

    return true;
}

template <typename T>
bool TernaryScheme<T>::isValidMerge(int p, const TernaryScheme<T> &scheme, int maxN, int maxRank) const {
    int j = (p + 1) % 3;
    int k = (p + 2) % 3;

    int maxElements = sizeof(T) * 8;
    int n = dimension[p] + scheme.dimension[p];

    return n <= maxN && n * dimension[j] <= maxElements && n * dimension[k] <= maxElements && dimension[j] == scheme.dimension[j] && dimension[k] == scheme.dimension[k] && rank + scheme.rank <= maxRank;
}

template <typename T>
bool TernaryScheme<T>::fixSigns() {
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
bool TernaryScheme<T>::validateDimensions() const {
    int maxSize = sizeof(T) * 8;

    for (int i = 0; i < 3; i++) {
        if (dimension[i] < 1 || dimension[i] > maxSize) {
            std::cout << "Invalid dimension \"" << dimension[i] << "\". Possible dimensions are 1 .. " << maxSize << std::endl;
            return false;
        }

        if (elements[i] < 1 || elements[i] > maxSize) {
            std::cout << "Invalid matrix elements count \"" << elements[i] << "\". Possible counts are 1 .. " << maxSize << std::endl;
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
bool TernaryScheme<T>::validateEquation(int i, int j, int k) const {
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
void TernaryScheme<T>::saveMatrix(std::ofstream &f, std::string name, const std::vector<TernaryVector<T>> &vectors) const {
    f << "    \"" << name << "\": [" << std::endl;

    for (size_t index = 0; index < vectors.size(); index++)
        f << "        [" << vectors[index] << "]" << (index < vectors.size() - 1 ? "," : "") << std::endl;

    f << "    ]";
}
