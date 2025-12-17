#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cassert>
#include <algorithm>

#include "../entities/flip_set.h"

template <typename T>
class BinaryScheme {
    int dimension[3];
    int elements[3];
    int rank;
    std::vector<T> uvw[3];
    FlipSet flips[3];

    std::uniform_int_distribution<int> boolDistribution;
    std::uniform_int_distribution<int> ijkDistribution;
public:
    BinaryScheme();
    BinaryScheme(const BinaryScheme<T> &scheme);

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
    bool tryReduce();

    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void copy(const BinaryScheme &scheme);

    bool validate() const;
private:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const T &u, const T &v, const T &w);

    void flip(int i, int j, int k, int index1, int index2);
    void plus(int i, int j, int k, int index1, int index2, int variant);
    void split(int i, int j, int k, int index1, int index2);
    void reduce(int i, int index1, int index2);
    bool checkFlipReduce(int j, int k, int index1, int index2);

    bool validateDimensions() const;
    bool validateEquation(int i, int j, int k) const;
    void saveMatrix(std::ofstream &f, std::string name, const std::vector<T> &vectors, int size) const;
};

template <typename T>
BinaryScheme<T>::BinaryScheme() : boolDistribution(0, 1), ijkDistribution(0, 2) {

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
bool BinaryScheme<T>::read(const std::string &path) {
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
            T vector = 0;
            for (int j = 0; j < elements[i]; j++) {
                int value;
                f >> value;
                vector |= T(abs(value) % 2) << j;
            }

            uvw[i].push_back(vector);
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
int BinaryScheme<T>::getRank() const {
    return rank;
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
int BinaryScheme<T>::getDimension(int index) const {
    return dimension[index];
}

template <typename T>
std::string BinaryScheme<T>::getRing() const {
    return "Z2";
}

template <typename T>
int BinaryScheme<T>::getAvailableFlips() const {
    return flips[0].size() + flips[1].size() + flips[2].size();
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
bool BinaryScheme<T>::validate() const {
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                if (!validateEquation(i, j, k))
                    return false;

    return true;
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
                f << ((uvw[i][index] >> j) & 1) << " ";

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
            f << (i > 0 ? ", " : "") << ((vectors[index] >> i) & 1);

        f << "]" << (index < vectors.size() - 1 ? "," : "") << std::endl;
    }

    f << "    ]";
}
