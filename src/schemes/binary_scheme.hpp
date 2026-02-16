#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cassert>
#include <algorithm>

#include "../algebra/binary_matrix.h"
#include "../algebra/binary_solver.h"
#include "../lift/binary_lifter.h"
#include "fractional_scheme.h"
#include "base_scheme.h"

template <typename T>
class BinaryScheme : public BaseScheme {
protected:
    std::vector<T> uvw[3];
public:
    BinaryScheme();
    BinaryScheme(const BinaryScheme<T> &scheme);

    bool initializeNaive(int n1, int n2, int n3);
    bool read(const std::string &path, bool checkCorrectness);
    bool read(std::istream &is, bool checkCorrectness);

    int getComplexity() const;
    std::string getRing() const;
    std::string getHash() const;

    bool tryFlip(std::mt19937 &generator);
    bool tryPlus(std::mt19937 &generator);
    bool trySplit(std::mt19937 &generator);
    bool tryExpand(std::mt19937 &generator);
    bool trySandwiching(std::mt19937 &generator);
    bool tryReduce();

    bool tryProject(std::mt19937 &generator, int minN);
    bool tryExtend(std::mt19937 &generator, int maxN, int maxRank);
    bool tryMerge(const BinaryScheme<T> &scheme, std::mt19937 &generator, int maxN, int maxRank);
    bool tryProduct(const BinaryScheme<T> &scheme, int maxN, int maxRank);

    void swapSizes(std::mt19937 &generator);
    void swapSizes(int p1, int p2);
    void merge(const BinaryScheme<T> &scheme, int p);
    void project(int p, int q);
    void extend(int p);
    void product(const BinaryScheme<T> &scheme);

    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void copy(const BinaryScheme &scheme);

    bool validate() const;
    bool reconstruct(FractionalScheme &scheme) const;

    BinaryLifter toLift() const;
protected:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const T &u, const T &v, const T &w);
    void excludeColumn(int matrix, int column);
    void excludeRow(int matrix, int row);
    void addColumn(int matrix);
    void addRow(int matrix);

    void flip(int i, int j, int k, int index1, int index2);
    void plus(int i, int j, int k, int index1, int index2, int variant);
    void split(int i, int j, int k, int index1, int index2);
    void reduce(int i, int index1, int index2);
    bool checkFlipReduce(int j, int k, int index1, int index2);

    bool isValidProject(int p, int minN) const;
    bool isValidExtension(int p, int maxN, int maxRank) const;
    bool isValidProduct(const BinaryScheme<T> &scheme, int maxN, int maxRank) const;
    bool isValidMerge(int p, const BinaryScheme<T> &scheme, int maxN, int maxRank) const;

    bool validateDimensions() const;
    bool validateEquation(int i, int j, int k) const;
    void saveMatrix(std::ofstream &f, std::string name, const std::vector<T> &vectors, int size) const;

    BinarySolver getJakobian() const;
};

template <typename T>
BinaryScheme<T>::BinaryScheme() {

}

template <typename T>
BinaryScheme<T>::BinaryScheme(const BinaryScheme &scheme) {
    copy(scheme);
}

template <typename T>
bool BinaryScheme<T>::initializeNaive(int n1, int n2, int n3) {
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
                uvw[0].push_back(T(1) << (i * n2 + k));
                uvw[1].push_back(T(1) << (k * n3 + j));
                uvw[2].push_back(T(1) << (j * n1 + i));
            }
        }
    }

    initFlips();
    return true;
}

template <typename T>
bool BinaryScheme<T>::read(const std::string &path, bool checkCorrectness) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    bool valid = read(f, checkCorrectness);
    f.close();

    if (!valid) {
        std::cout << "Invalid scheme in the file \"" << path << "\"" << std::endl;
        return false;
    }

    return true;
}

template <typename T>
bool BinaryScheme<T>::read(std::istream &is, bool checkCorrectness) {
    is >> dimension[0] >> dimension[1] >> dimension[2] >> rank;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    if (!validateDimensions())
        return false;

    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++) {
            T vector = 0;
            for (int j = 0; j < elements[i]; j++) {
                int value;
                is >> value;
                vector |= T(abs(value) % 2) << j;
            }

            uvw[i].push_back(vector);
        }
    }

    if (checkCorrectness && !validate())
        return false;

    initFlips();
    return true;
}

template <typename T>
int BinaryScheme<T>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += __builtin_popcountll(uvw[i][index]);

    return count - 2 * rank - elements[2];
}

template <typename T>
std::string BinaryScheme<T>::getRing() const {
    return "Z2";
}

template <typename T>
std::string BinaryScheme<T>::getHash() const {
    std::vector<std::string> lines;

    for (int index = 0; index < rank; index++) {
        std::stringstream ss;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < elements[i]; j++)
                ss << (int)((uvw[i][index] >> j) & 1);

        lines.push_back(ss.str());
    }

    std::sort(lines.begin(), lines.end());
    std::stringstream hash;
    for (int index = 0; index < rank; index++)
        hash << lines[index];

    return hash.str();
}

template <typename T>
bool BinaryScheme<T>::tryFlip(std::mt19937 &generator) {
    size_t size = flips[0].size() + flips[1].size() + flips[2].size();

    if (!size)
        return false;

    size_t index = generator() % size;

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

    flip(i, j, k, index1, index2);
    return true;
}

template <typename T>
bool BinaryScheme<T>::tryPlus(std::mt19937 &generator) {
    std::uniform_int_distribution<int> distribution(0, rank - 1);
    int index1 = distribution(generator);
    int index2 = distribution(generator);

    while (index1 == index2 || uvw[0][index1] == uvw[0][index2] || uvw[1][index1] == uvw[1][index2] || uvw[2][index1] == uvw[2][index2]) {
        index1 = distribution(generator);
        index2 = distribution(generator);
    }

    int permutation[3] = {0, 1, 2};
    std::shuffle(permutation, permutation + 3, generator);

    plus(permutation[0], permutation[1], permutation[2], index1, index2, ijkDistribution(generator));
    return true;
}

template <typename T>
bool BinaryScheme<T>::trySplit(std::mt19937 &generator) {
    std::uniform_int_distribution<int> distribution(0, rank - 1);
    int index1, index2, i;

    do {
        index1 = distribution(generator);
        index2 = distribution(generator);
        i = ijkDistribution(generator);
    } while (index1 == index2 || uvw[i][index1] == uvw[i][index2]);

    split(i, (i + 1) % 3, (i + 2) % 3, index1, index2);
    return true;
}

template <typename T>
bool BinaryScheme<T>::tryExpand(std::mt19937 &generator) {
    if (rank >= dimension[0] * dimension[1] * dimension[2])
        return false;

    if (boolDistribution(generator))
        return tryPlus(generator);

    return trySplit(generator);
}

template <typename T>
bool BinaryScheme<T>::trySandwiching(std::mt19937 &generator) {
    BinaryMatrix u(dimension[0], dimension[0]);
    BinaryMatrix v(dimension[1], dimension[1]);
    BinaryMatrix w(dimension[2], dimension[2]);

    BinaryMatrix u1(dimension[0], dimension[0]);
    BinaryMatrix v1(dimension[1], dimension[1]);
    BinaryMatrix w1(dimension[2], dimension[2]);

    u.randomInvertible(u1, generator);
    v.randomInvertible(v1, generator);
    w.randomInvertible(w1, generator);

    for (int index = 0; index < rank; index++) {
        BinaryMatrix mu(dimension[0], dimension[1]);
        BinaryMatrix mv(dimension[1], dimension[2]);
        BinaryMatrix mw(dimension[2], dimension[0]);

        for (int i = 0; i < elements[0]; i++)
            mu[i] = (uvw[0][index] >> i) & 1;

        for (int i = 0; i < elements[1]; i++)
            mv[i] = (uvw[1][index] >> i) & 1;

        for (int i = 0; i < elements[2]; i++)
            mw[i] = (uvw[2][index] >> i) & 1;

        mu.sandwich(u, v1);
        mv.sandwich(v, w1);
        mw.sandwich(w, u1);

        uvw[0][index] = 0;
        for (int i = 0; i < elements[0]; i++)
            uvw[0][index] |= T(mu[i]) << i;

        uvw[1][index] = 0;
        for (int i = 0; i < elements[1]; i++)
            uvw[1][index] |= T(mv[i]) << i;

        uvw[2][index] = 0;
        for (int i = 0; i < elements[2]; i++)
            uvw[2][index] |= T(mw[i]) << i;
    }

    initFlips();
    return true;
}

template <typename T>
bool BinaryScheme<T>::tryReduce() {
    for (size_t i = 0; i < flips[0].size(); i++) {
        int index1 = flips[0].index1(i);
        int index2 = flips[0].index2(i);

        if (uvw[1][index1] == uvw[1][index2]) {
            reduce(2, index1, index2);
            return true;
        }

        if (uvw[2][index1] == uvw[2][index2]) {
            reduce(1, index1, index2);
            return true;
        }
    }

    for (size_t i = 0; i < flips[1].size(); i++) {
        int index1 = flips[1].index1(i);
        int index2 = flips[1].index2(i);

        if (uvw[2][index1] == uvw[2][index2]) {
            reduce(0, index1, index2);
            return true;
        }
    }

    return false;
}

template <typename T>
bool BinaryScheme<T>::tryProject(std::mt19937 &generator, int minN) {
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
bool BinaryScheme<T>::tryExtend(std::mt19937 &generator, int maxN, int maxRank) {
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
bool BinaryScheme<T>::tryMerge(const BinaryScheme<T> &scheme, std::mt19937 &generator, int maxN, int maxRank) {
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
bool BinaryScheme<T>::tryProduct(const BinaryScheme<T> &scheme, int maxN, int maxRank) {
    if (!isValidProduct(scheme, maxN, maxRank))
        return false;

    product(scheme);
    return true;
}

template <typename T>
void BinaryScheme<T>::swapSizes(std::mt19937 &generator) {
    int p1, p2;

    do {
        p1 = ijkDistribution(generator);
        p2 = ijkDistribution(generator);
    } while (p1 == p2);

    swapSizes(p1, p2);
}

template <typename T>
void BinaryScheme<T>::swapSizes(int p1, int p2) {
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
        T u = 0;;
        T v = 0;
        T w = 0;

        for (int i = 0; i < dimensionNew[0]; i++)
            for (int j = 0; j < dimensionNew[1]; j++)
                u |= T((uvw[indices[0]][index] >> (j * dimensionNew[0] + i)) & 1) << (i * dimensionNew[1] + j);

        for (int i = 0; i < dimensionNew[1]; i++)
            for (int j = 0; j < dimensionNew[2]; j++)
                v |= T((uvw[indices[1]][index] >> (j * dimensionNew[1] + i)) & 1) << (i * dimensionNew[2] + j);

        for (int i = 0; i < dimensionNew[2]; i++)
            for (int j = 0; j < dimensionNew[0]; j++)
                w |= T((uvw[indices[2]][index] >> (j * dimensionNew[2] + i)) & 1) << (i * dimensionNew[0] + j);

        uvw[0][index] = u;
        uvw[1][index] = v;
        uvw[2][index] = w;
    }

    for (int i = 0; i < 3; i++) {
        dimension[i] = dimensionNew[i];
        elements[i] = dimensionNew[i] * dimensionNew[(i + 1) % 3];
    }

    initFlips();
}

template <typename T>
void BinaryScheme<T>::merge(const BinaryScheme<T> &scheme, int p) {
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
        T u = 0;
        T v = 0;
        T w = 0;

        for (int i = 0; i < dimension[0]; i++)
            for (int j = 0; j < dimension[1]; j++)
                u |= T((uvw[0][index] >> (i * dimension[1] + j)) & 1) << (i * dimensionNew[1] + j);

        for (int i = 0; i < dimension[1]; i++)
            for (int j = 0; j < dimension[2]; j++)
                v |= T((uvw[1][index] >> (i * dimension[2] + j)) & 1) << (i * dimensionNew[2] + j);

        for (int i = 0; i < dimension[2]; i++)
            for (int j = 0; j < dimension[0]; j++)
                w |= T((uvw[2][index] >> (i * dimension[0] + j)) & 1) << (i * dimensionNew[0] + j);

        uvw[0][index] = u;
        uvw[1][index] = v;
        uvw[2][index] = w;
    }

    for (int index = 0; index < scheme.rank; index++) {
        T u = 0;
        T v = 0;
        T w = 0;

        for (int i = 0; i < scheme.dimension[0]; i++)
            for (int j = 0; j < scheme.dimension[1]; j++)
                u |= T((scheme.uvw[0][index] >> (i * scheme.dimension[1] + j)) & 1) << ((i + d[0]) * dimensionNew[1] + j + d[1]);

        for (int i = 0; i < scheme.dimension[1]; i++)
            for (int j = 0; j < scheme.dimension[2]; j++)
                v |= T((scheme.uvw[1][index] >> (i * scheme.dimension[2] + j)) & 1) << ((i + d[1]) * dimensionNew[2] + j + d[2]);

        for (int i = 0; i < scheme.dimension[2]; i++)
            for (int j = 0; j < scheme.dimension[0]; j++)
                w |= T((scheme.uvw[2][index] >> (i * scheme.dimension[0] + j)) & 1) << ((i + d[2]) * dimensionNew[0] + j + d[0]);

        addTriplet(0, 1, 2, u, v, w);
    }

    for (int i = 0; i < 3; i++) {
        dimension[i] = dimensionNew[i];
        elements[i] = elementsNew[i];
    }

    initFlips();
}

template <typename T>
void BinaryScheme<T>::project(int p, int q) {
    excludeRow(p, q);
    excludeColumn((p + 2) % 3, q);
    dimension[p]--;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    removeZeroes();
    initFlips();
}

template <typename T>
void BinaryScheme<T>::extend(int p) {
    addRow(p);
    addColumn((p + 2) % 3);

    if (p == 0) {
        for (int i = 0; i < dimension[2]; i++) {
            for (int j = 0; j < dimension[1]; j++) {
                T u = T(1) << (dimension[0] * dimension[1] + j);
                T v = T(1) << (j * dimension[2] + i);
                T w = T(1) << (i * (dimension[0] + 1) + dimension[0]);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }
    else if (p == 1) {
        for (int i = 0; i < dimension[0]; i++) {
            for (int j = 0; j < dimension[2]; j++) {
                T u = T(1) << (i * (dimension[1] + 1) + dimension[1]);
                T v = T(1) << (dimension[1] * dimension[2] + j);
                T w = T(1) << (j * dimension[0] + i);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }
    else {
        for (int i = 0; i < dimension[0]; i++) {
            for (int j = 0; j < dimension[1]; j++) {
                T u = T(1) << (i * dimension[1] + j);
                T v = T(1) << (j * (dimension[2] + 1) + dimension[2]);
                T w = T(1) << (dimension[2] * dimension[0] + i);
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
void BinaryScheme<T>::product(const BinaryScheme<T> &scheme2) {
    BinaryScheme<T> scheme1(*this);

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

                T vector = 0;

                for (int i = 0; i < scheme1.elements[p]; i++) {
                    for (int j = 0; j < scheme2.elements[p]; j++) {
                        int row1 = i / scheme1.dimension[p1];
                        int col1 = i % scheme1.dimension[p1];
                        T value1 = (scheme1.uvw[p][index1] >> i) & 1;

                        int row2 = j / scheme2.dimension[p1];
                        int col2 = j % scheme2.dimension[p1];
                        T value2 = (scheme2.uvw[p][index2] >> j) & 1;

                        int row = row1 * scheme2.dimension[p] + row2;
                        int col = col1 * scheme2.dimension[p1] + col2;

                        vector |= T(value1 * value2) << (row * dimension[p1] + col);
                    }
                }

                uvw[p].push_back(vector);
            }
        }
    }

    initFlips();
}

template <typename T>
void BinaryScheme<T>::saveJson(const std::string &path) const {
    std::ofstream f(path);

    f << "{" << std::endl;
    f << "    \"n\": [" << dimension[0] << ", " << dimension[1] << ", " << dimension[2] << "]," << std::endl;
    f << "    \"m\": " << rank << "," << std::endl;
    f << "    \"z2\": true," << std::endl;
    f << "    \"complexity\": " << getComplexity() << "," << std::endl;

    saveMatrix(f, "u", uvw[0], elements[0]);
    f << "," << std::endl;
    saveMatrix(f, "v", uvw[1], elements[1]);
    f << "," << std::endl;
    saveMatrix(f, "w", uvw[2], elements[2]);
    f << std::endl;
    f << "}" << std::endl;

    f.close();
}

template <typename T>
void BinaryScheme<T>::saveTxt(const std::string &path) const {
    std::ofstream f(path);

    f << dimension[0] << " " << dimension[1] << " " << dimension[2] << " " << rank << std::endl;
    
    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++)
            for (int j = 0; j < elements[i]; j++)
                f << (int)((uvw[i][index] >> j) & 1) << " ";

        f << std::endl;
    }

    f.close();
}

template <typename T>
void BinaryScheme<T>::copy(const BinaryScheme &scheme) {
    rank = scheme.rank;

    for (int i = 0; i < 3; i++) {
        dimension[i] = scheme.dimension[i];
        elements[i] = scheme.elements[i];
        uvw[i].clear();

        for (int index = 0; index < rank; index++)
            uvw[i].push_back(scheme.uvw[i][index]);
    }

    initFlips();
}

template <typename T>
bool BinaryScheme<T>::validate() const {
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                if (!validateEquation(i, j, k))
                    return false;

    return true;
}

template <typename T>
bool BinaryScheme<T>::reconstruct(FractionalScheme &scheme) const {
    return false;
}

template <typename T>
BinaryLifter BinaryScheme<T>::toLift() const {
    std::vector<uint64_t> u(rank * elements[0]);
    std::vector<uint64_t> v(rank * elements[1]);
    std::vector<uint64_t> w(rank * elements[2]);

    for (int index = 0; index < rank; index++) {
        for (int i = 0; i < elements[0]; i++)
            u[index * elements[0] + i] = (uvw[0][index] >> i) & 1;

        for (int i = 0; i < elements[1]; i++)
            v[index * elements[1] + i] = (uvw[1][index] >> i) & 1;

        for (int i = 0; i < elements[2]; i++)
            w[index * elements[2] + i] = (uvw[2][index] >> i) & 1;
    }

    return BinaryLifter(dimension[0], dimension[1], dimension[2], rank, u, v, w, getJakobian());
}

template <typename T>
void BinaryScheme<T>::initFlips() {
    for (int i = 0; i < 3; i++) {
        flips[i].clear();

        for (int index1 = 0; index1 < rank; index1++)
            for (int index2 = index1 + 1; index2 < rank; index2++)
                if (uvw[i][index1] == uvw[i][index2])
                    flips[i].add(index1, index2);
    }
}

template <typename T>
void BinaryScheme<T>::removeZeroes() {
    for (int index = 0; index < rank; index++)
        if (!uvw[0][index] || !uvw[1][index] || !uvw[2][index])
            removeAt(index--);
}

template <typename T>
void BinaryScheme<T>::removeAt(int index) {
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
void BinaryScheme<T>::addTriplet(int i, int j, int k, const T &u, const T &v, const T &w) {
    uvw[i].push_back(u);
    uvw[j].push_back(v);
    uvw[k].push_back(w);
    rank++;
}

template <typename T>
void BinaryScheme<T>::excludeColumn(int matrix, int column) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];
    std::vector<int> oldColumns;

    for (int j = 0; j < n2; j++)
        if (j != column)
            oldColumns.push_back(j);

    for (int index = 0; index < rank; index++) {
        T vector = 0;

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2 - 1; j++)
                vector |= T((uvw[matrix][index] >> (i * n2 + oldColumns[j])) & 1) << (i * (n2 - 1) + j);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void BinaryScheme<T>::excludeRow(int matrix, int row) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];
    std::vector<int> oldRows;

    for (int i = 0; i < n1; i++)
        if (i != row)
            oldRows.push_back(i);

    for (int index = 0; index < rank; index++) {
        T vector = 0;

        for (int i = 0; i < n1 - 1; i++)
            for (int j = 0; j < n2; j++)
                vector |= T((uvw[matrix][index] >> (oldRows[i] * n2 + j)) & 1) << (i * n2 + j);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void BinaryScheme<T>::addColumn(int matrix) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];

    for (int index = 0; index < rank; index++) {
        T vector = 0;

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vector |= T((uvw[matrix][index] >> (i * n2 + j)) & 1) << (i * (n2 + 1) + j);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void BinaryScheme<T>::addRow(int matrix) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];

    for (int index = 0; index < rank; index++) {
        T vector = 0;

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vector |= T((uvw[matrix][index] >> (i * n2 + j)) & 1) << (i * n2 + j);

        uvw[matrix][index] = vector;
    }
}

template <typename T>
void BinaryScheme<T>::flip(int i, int j, int k, int index1, int index2) {
    uvw[j][index1] ^= uvw[j][index2];
    uvw[k][index2] ^= uvw[k][index1];

    flips[j].remove(index1);
    flips[k].remove(index2);

    if (!uvw[j][index1] || !uvw[k][index2]) {
        removeZeroes();
        initFlips();
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
void BinaryScheme<T>::plus(int i, int j, int k, int index1, int index2, int variant) {
    const T a1 = uvw[i][index1];
    const T b1 = uvw[j][index1];
    const T c1 = uvw[k][index1];

    const T a2 = uvw[i][index2];
    const T b2 = uvw[j][index2];
    const T c2 = uvw[k][index2];

    const T a = a1 ^ a2;
    const T b = b1 ^ b2;
    const T c = c1 ^ c2;

    if (variant == 0) {
        uvw[j][index1] = b;
        uvw[i][index2] = a;
        addTriplet(i, j, k, a1, b2, c);
    }
    else if (variant == 1) {
        uvw[k][index1] = c;
        uvw[j][index2] = b;
        addTriplet(i, j, k, a, b1, c2);
    }
    else {
        uvw[i][index1] = a;
        uvw[k][index2] = c;
        addTriplet(i, j, k, a2, b, c1);
    }

    if (!a || !b || !c)
        removeZeroes();

    initFlips();
}

template <typename T>
void BinaryScheme<T>::split(int i, int j, int k, int index1, int index2) {
    const T u = uvw[i][index1] ^ uvw[i][index2];
    const T v = uvw[j][index1];
    const T w = uvw[k][index1];

    addTriplet(i, j, k, u, v, w);
    uvw[i][index1] = uvw[i][index2];

    removeZeroes();
    initFlips();
}

template <typename T>
void BinaryScheme<T>::reduce(int i, int index1, int index2) {
    uvw[i][index1] ^= uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <typename T>
bool BinaryScheme<T>::checkFlipReduce(int i, int j, int index1, int index2) {
    if (uvw[i][index1] == uvw[i][index2]) {
        reduce(j, index1, index2);
        return true;
    }

    if (uvw[j][index1] == uvw[j][index2]) {
        reduce(i, index1, index2);
        return true;
    }

    return false;
}


template <typename T>
bool BinaryScheme<T>::isValidProject(int p, int minN) const {
    return dimension[p] > minN && dimension[(p + 1) % 3] >= minN && dimension[(p + 2) % 3] >= minN;
}

template <typename T>
bool BinaryScheme<T>::isValidExtension(int p, int maxN, int maxRank) const {
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
bool BinaryScheme<T>::isValidProduct(const BinaryScheme<T> &scheme, int maxN, int maxRank) const {
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
bool BinaryScheme<T>::isValidMerge(int p, const BinaryScheme<T> &scheme, int maxN, int maxRank) const {
    int j = (p + 1) % 3;
    int k = (p + 2) % 3;

    int maxElements = sizeof(T) * 8;
    int n = dimension[p] + scheme.dimension[p];

    return n <= maxN && n * dimension[j] <= maxElements && n * dimension[k] <= maxElements && dimension[j] == scheme.dimension[j] && dimension[k] == scheme.dimension[k] && rank + scheme.rank <= maxRank;
}

template <typename T>
bool BinaryScheme<T>::validateDimensions() const {
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
bool BinaryScheme<T>::validateEquation(int i, int j, int k) const {
    int i1 = i / dimension[1];
    int i2 = i % dimension[1];
    int j1 = j / dimension[2];
    int j2 = j % dimension[2];
    int k1 = k / dimension[0];
    int k2 = k % dimension[0];

    bool target = (i2 == j1) && (i1 == k2) && (j2 == k1);
    bool equation = false;

    for (int index = 0; index < rank; index++)
        equation ^= ((uvw[0][index] >> i) & 1) && ((uvw[1][index] >> j) & 1) && ((uvw[2][index] >> k) & 1);

    return equation == target;
}

template <typename T>
void BinaryScheme<T>::saveMatrix(std::ofstream &f, std::string name, const std::vector<T> &vectors, int size) const {
    f << "    \"" << name << "\": [" << std::endl;

    for (size_t index = 0; index < vectors.size(); index++) {
        f << "        [";

        for (int i = 0; i < size; i++)
            f << (i > 0 ? ", " : "") << (int)((vectors[index] >> i) & 1);

        f << "]" << (index < vectors.size() - 1 ? "," : "") << std::endl;
    }

    f << "    ]";
}

template <typename T>
BinarySolver BinaryScheme<T>::getJakobian() const {
    int rows = elements[0] * elements[1] * elements[2];
    int columns = rank * (elements[0] + elements[1] + elements[2]);
    BinarySolver jakobian(rows, columns);

    int vOffset = elements[0] * rank;
    int wOffset = (elements[0] + elements[1]) * rank;

    for (int i = 0; i < elements[0]; i++) {
        for (int j = 0; j < elements[1]; j++) {
            for (int k = 0; k < elements[2]; k++) {
                int row = (i * elements[1] + j) * elements[2] + k;

                for (int index = 0; index < rank; index++) {
                    uint8_t u = (uvw[0][index] >> i) & 1;
                    uint8_t v = (uvw[1][index] >> j) & 1;
                    uint8_t w = (uvw[2][index] >> k) & 1;

                    jakobian.set(row, i * rank + index, v & w);
                    jakobian.set(row, vOffset + j * rank + index, u & w);
                    jakobian.set(row, wOffset + k * rank + index, u & v);
                }
            }
        }
    }

    return jakobian;
}
