#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cassert>
#include <algorithm>

#include "../entities/mod_vector.h"
#include "../entities/ranks.h"
#include "../entities/invariants_builder.h"
#include "../algebra/matrix.h"
#include "../algebra/mod_matrix.h"
#include "../lift/mod_lifter.h"
#include "fractional_scheme.h"
#include "base_scheme.h"

template <int64_t M>
class ModScheme : public BaseScheme {
protected:
    std::vector<ModVector> uvw[3];
private:
    FlipSet flipsNeg[3];
public:
    ModScheme();
    ModScheme(const ModScheme<M> &scheme);

    bool initializeNaive(int n1, int n2, int n3);
    bool read(const std::string &path, bool checkCorrectness);
    bool read(std::istream &is, bool checkCorrectness);

    int getComplexity() const;
    std::string getRing() const;
    std::string getHash() const;

    std::vector<Ranks> getRanks() const;
    std::string getTypeInvariant() const;
    std::string getMinTypeInvariant() const;

    std::vector<std::vector<int>> getExpressionsU() const;
    std::vector<std::vector<int>> getExpressionsV() const;
    std::vector<std::vector<int>> getExpressionsW() const;

    bool tryFlip(std::mt19937 &generator);
    bool tryPlus(std::mt19937 &generator);
    bool trySplit(std::mt19937 &generator);
    bool tryExpand(std::mt19937 &generator);
    bool trySandwiching(std::mt19937 &generator);
    bool tryReduce();
    bool check2Reduce() const;

    bool tryProject(std::mt19937 &generator, int minN);
    bool tryExtend(std::mt19937 &generator, int maxN, int maxRank);
    bool tryMerge(const ModScheme<M> &scheme, std::mt19937 &generator, int maxN, int maxRank);
    bool tryProduct(const ModScheme<M> &scheme, int maxN, int maxRank);

    bool isValidProject(int p, int minN) const;
    bool isValidExtension(int p, int maxN, int maxRank) const;
    bool isValidProduct(const ModScheme<M> &scheme, int maxN, int maxRank) const;
    bool isValidMerge(int p, const ModScheme<M> &scheme, int maxN, int maxRank) const;

    void fixSizes();
    void setSizes(int n1, int n2, int n3);
    void swapSizes(std::mt19937 &generator);
    void swapSizes(int p1, int p2);
    void merge(const ModScheme<M> &scheme, int p);
    void project(int p, int q);
    void extend(int p);
    void product(const ModScheme<M> &scheme);

    void saveJson(const std::string &path) const;
    void saveTxt(const std::string &path) const;
    void copy(const ModScheme &scheme);

    bool validate() const;
    bool validateParallel() const;
    bool reconstruct(FractionalScheme &scheme) const;

    ModLifter toLift() const;
    bool lift(FractionalScheme& lifted, int steps) const;
    bool canLift(int steps) const;
private:
    void initFlips();
    void removeZeroes();
    void removeAt(int index);
    void addTriplet(int i, int j, int k, const ModVector &u, const ModVector &v, const ModVector &w);
    void excludeColumn(int matrix, int column);
    void excludeRow(int matrix, int row);
    void addColumn(int matrix);
    void addRow(int matrix);

    void flip(int i, int j, int k, int index1, int index2, bool pos);
    void plus(int i, int j, int k, int index1, int index2, int variant);
    void split(int i, int j, int k, int index1, int index2);
    bool reduce(int i, int index1, int index2, int sign);
    void reduceAdd(int i, int index1, int index2);
    void reduceSub(int i, int index1, int index2);
    bool checkFlipReduce(int j, int k, int index1, int index2, int sign);


    bool validateDimensions() const;
    bool validateEquation(int i, int j, int k) const;
    void normalize();
    void saveMatrix(std::ofstream &f, std::string name, const std::vector<ModVector> &vectors) const;

    ModMatrix getJakobian() const;
};

template <int64_t M>
ModScheme<M>::ModScheme() {

}

template <int64_t M>
ModScheme<M>::ModScheme(const ModScheme &scheme) {
    copy(scheme);
}

template <int64_t M>
bool ModScheme<M>::initializeNaive(int n1, int n2, int n3) {
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
                uvw[0].emplace_back(ModVector(n1 * n2, i * n2 + k, M));
                uvw[1].emplace_back(ModVector(n2 * n3, k * n3 + j, M));
                uvw[2].emplace_back(ModVector(n3 * n1, j * n1 + i, M));
            }
        }
    }

    initFlips();
    return true;
}

template <int64_t M>
bool ModScheme<M>::read(const std::string &path, bool checkCorrectness) {
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

template <int64_t M>
bool ModScheme<M>::read(std::istream &is, bool checkCorrectness) {
    is >> dimension[0] >> dimension[1] >> dimension[2] >> rank;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    if (!validateDimensions())
        return false;

    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++) {
            ModVector vector(elements[i], M);
            is >> vector;
            uvw[i].emplace_back(vector);
        }
    }

    if (checkCorrectness && !validate())
        return false;

    initFlips();
    return true;
}

template <int64_t M>
int ModScheme<M>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += uvw[i][index].nonZeroCount();

    return count - 2 * rank - elements[2];
}

template <int64_t M>
std::string ModScheme<M>::getRing() const {
    return "Z" + std::to_string(M);
}

template <int64_t M>
std::string ModScheme<M>::getHash() const {
    std::vector<std::string> lines;

    for (int index = 0; index < rank; index++) {
        std::stringstream ss;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < elements[i]; j++)
                ss << uvw[i][index][j];

        lines.push_back(ss.str());
    }

    std::sort(lines.begin(), lines.end());
    std::stringstream hash;
    for (int index = 0; index < rank; index++)
        hash << lines[index];

    return hash.str();
}

template <int64_t M>
std::vector<Ranks> ModScheme<M>::getRanks() const {
    std::vector<Ranks> ranks;

    for (int index = 0; index < rank; index++) {
        ModMatrix u(dimension[0], dimension[1], M);
        ModMatrix v(dimension[1], dimension[2], M);
        ModMatrix w(dimension[2], dimension[0], M);

        for (int i = 0; i < elements[0]; i++)
            u[i] = uvw[0][index][i];

        for (int i = 0; i < elements[1]; i++)
            v[i] = uvw[1][index][i];

        for (int i = 0; i < elements[2]; i++)
            w[i] = uvw[2][index][i];

        ranks.push_back({u.rank(), v.rank(), w.rank()});
    }

    return ranks;
}

template <int64_t M>
std::string ModScheme<M>::getTypeInvariant() const {
    std::vector<Ranks> ranks = getRanks();
    InvariantsBuilder invariants(ranks);
    return invariants.getType();
}

template <int64_t M>
std::string ModScheme<M>::getMinTypeInvariant() const {
    std::vector<Ranks> ranks = getRanks();
    InvariantsBuilder invariants(ranks);
    return invariants.getMinType();
}

template <int64_t M>
std::vector<std::vector<int>> ModScheme<M>::getExpressionsU() const {
    std::vector<std::vector<int>> expressions(rank, std::vector<int>(elements[0], 0));

    for (int index = 0; index < rank; index++)
        for (int i = 0; i < elements[0]; i++)
            expressions[index][i] = uvw[0][index][i];

    return expressions;
}

template <int64_t M>
std::vector<std::vector<int>> ModScheme<M>::getExpressionsV() const {
    std::vector<std::vector<int>> expressions(rank, std::vector<int>(elements[1], 0));

    for (int index = 0; index < rank; index++)
        for (int i = 0; i < elements[1]; i++)
            expressions[index][i] = uvw[1][index][i];

    return expressions;
}

template <int64_t M>
std::vector<std::vector<int>> ModScheme<M>::getExpressionsW() const {
    std::vector<std::vector<int>> expressions(elements[2], std::vector<int>(rank, 0));

    for (int i = 0; i < elements[2]; i++)
        for (int index = 0; index < rank; index++)
            expressions[i][index] = uvw[2][index][i];

    return expressions;
}

template <int64_t M>
bool ModScheme<M>::tryFlip(std::mt19937 &generator) {
    size_t sizePos = flips[0].size() + flips[1].size() + flips[2].size();
    size_t sizeNeg = flipsNeg[0].size() + flipsNeg[1].size() + flipsNeg[2].size();
    size_t size = sizePos + sizeNeg;

    if (!size)
        return false;

    size_t index = generator() % size;
    bool pos = index < sizePos;

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
    else if (index < sizePos) {
        i = 2;
        j = 0;
        k = 1;
        index -= flips[0].size() + flips[1].size();
    }
    else if (index < sizePos + flipsNeg[0].size()) {
        i = 0;
        j = 1;
        k = 2;
        index -= sizePos;
    }
    else if (index < sizePos + flipsNeg[0].size() + flipsNeg[1].size()) {
        i = 1;
        j = 0;
        k = 2;
        index -= sizePos + flipsNeg[0].size();
    }
    else {
        i = 2;
        j = 0;
        k = 1;
        index -= sizePos + flipsNeg[0].size() + flipsNeg[1].size();
    }

    int index1 = pos ? flips[i].index1(index) : flipsNeg[i].index1(index);
    int index2 = pos ? flips[i].index2(index) : flipsNeg[i].index2(index);

    if (boolDistribution(generator))
        std::swap(j, k);

    if (boolDistribution(generator))
        std::swap(index1, index2);

    flip(i, j, k, index1, index2, pos);
    return true;
}

template <int64_t M>
bool ModScheme<M>::tryPlus(std::mt19937 &generator) {
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

template <int64_t M>
bool ModScheme<M>::trySplit(std::mt19937 &generator) {
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

template <int64_t M>
bool ModScheme<M>::tryExpand(std::mt19937 &generator) {
    if (rank >= dimension[0] * dimension[1] * dimension[2])
        return false;

    if (boolDistribution(generator))
        return tryPlus(generator);

    return trySplit(generator);
}

template <int64_t M>
bool ModScheme<M>::trySandwiching(std::mt19937 &generator) {
    Matrix u(dimension[0], dimension[0]);
    Matrix v(dimension[1], dimension[1]);
    Matrix w(dimension[2], dimension[2]);

    Matrix u1(dimension[0], dimension[0]);
    Matrix v1(dimension[1], dimension[1]);
    Matrix w1(dimension[2], dimension[2]);

    do {
        u.random(0, M - 1, 1, generator);
    } while (!u.invertible(u1));

    do {
        v.random(0, M - 1, 1, generator);
    } while (!v.invertible(v1));

    do {
        w.random(0, M - 1, 1, generator);
    } while (!w.invertible(w1));

    if (!u1.toRing(M) || !v1.toRing(M) || !w1.toRing(M))
        return false;

    for (int index = 0; index < rank; index++) {
        Matrix mu(dimension[0], dimension[1]);
        Matrix mv(dimension[1], dimension[2]);
        Matrix mw(dimension[2], dimension[0]);

        for (int i = 0; i < elements[0]; i++)
            mu[i] = uvw[0][index][i];

        for (int i = 0; i < elements[1]; i++)
            mv[i] = uvw[1][index][i];

        for (int i = 0; i < elements[2]; i++)
            mw[i] = uvw[2][index][i];

        mu.sandwich(u, v1);
        mv.sandwich(v, w1);
        mw.sandwich(w, u1);

        mu.toRing(M);
        mv.toRing(M);
        mw.toRing(M);

        for (int i = 0; i < elements[0]; i++)
            uvw[0][index].set(i, mu[i].numerator());

        for (int i = 0; i < elements[1]; i++)
            uvw[1][index].set(i, mv[i].numerator());

        for (int i = 0; i < elements[2]; i++)
            uvw[2][index].set(i, mw[i].numerator());
    }

    initFlips();
    return true;
}

template <int64_t M>
bool ModScheme<M>::tryReduce() {
    for (size_t i = 0; i < flips[0].size(); i++) {
        int index1 = flips[0].index1(i);
        int index2 = flips[0].index2(i);

        if (reduce(2, index1, index2, uvw[1][index1].compare(uvw[1][index2])))
            return true;

        if (reduce(1, index1, index2, uvw[2][index1].compare(uvw[2][index2])))
            return true;
    }

    for (size_t i = 0; i < flips[1].size(); i++) {
        int index1 = flips[1].index1(i);
        int index2 = flips[1].index2(i);

        if (reduce(0, index1, index2, uvw[2][index1].compare(uvw[2][index2])))
            return true;
    }

    for (size_t i = 0; i < flipsNeg[0].size(); i++) {
        int index1 = flipsNeg[0].index1(i);
        int index2 = flipsNeg[0].index2(i);

        if (reduce(2, index1, index2, -uvw[1][index1].compare(uvw[1][index2])))
            return true;

        if (reduce(1, index1, index2, -uvw[2][index1].compare(uvw[2][index2])))
            return true;
    }

    for (size_t i = 0; i < flipsNeg[1].size(); i++) {
        int index1 = flipsNeg[1].index1(i);
        int index2 = flipsNeg[1].index2(i);

        if (reduce(0, index1, index2, -uvw[2][index1].compare(uvw[2][index2])))
            return true;
    }

    return false;
}

template <int64_t M>
bool ModScheme<M>::check2Reduce() const {
    for (int p = 0; p < 3; p++) {
        int p1 = (p + 1) % 3;
        int p2 = (p + 2) % 3;

        ModMatrix matrix(rank, elements[p1] * elements[p2], M);

        for (int index = 0; index < rank; index++)
            for (int i = 0; i < elements[p1]; i++)
                for (int j = 0; j < elements[p2]; j++)
                    matrix(index, i * elements[p2] + j) = (uvw[p1][index][i] * uvw[p2][index][j] + M) % M;

        if (matrix.rank() < rank)
            return true;
    }

    return false;
}

template <int64_t M>
bool ModScheme<M>::tryProject(std::mt19937 &generator, int minN) {
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

template <int64_t M>
bool ModScheme<M>::tryExtend(std::mt19937 &generator, int maxN, int maxRank) {
    std::vector<int> indices;

    for (int i = 0; i < 3; i++)
        if (isValidExtension(i, maxN, maxRank))
            indices.push_back(i);

    if (!indices.size())
        return false;

    extend(indices[generator() % indices.size()]);
    return true;
}

template <int64_t M>
bool ModScheme<M>::tryMerge(const ModScheme<M> &scheme, std::mt19937 &generator, int maxN, int maxRank) {
    std::vector<int> indices;

    for (int i = 0; i < 3; i++)
        if (isValidMerge(i, scheme, maxN, maxRank))
            indices.push_back(i);

    if (!indices.size())
        return false;

    merge(scheme, indices[generator() % indices.size()]);
    return true;
}

template <int64_t M>
bool ModScheme<M>::tryProduct(const ModScheme<M> &scheme, int maxN, int maxRank) {
    if (!isValidProduct(scheme, maxN, maxRank))
        return false;

    product(scheme);
    return true;
}

template <int64_t M>
void ModScheme<M>::fixSizes() {
    if (dimension[0] > dimension[1])
        swapSizes(0, 1);

    if (dimension[1] > dimension[2])
        swapSizes(1, 2);

    if (dimension[0] > dimension[1])
        swapSizes(0, 1);
}

template <int64_t M>
void ModScheme<M>::setSizes(int n1, int n2, int n3) {
    if (dimension[0] == n1 && dimension[1] == n3 && dimension[2] == n2) {
        swapSizes(1, 2);
    }
    else if (dimension[0] == n2 && dimension[1] == n1 && dimension[2] == n3) {
        swapSizes(0, 1);
    }
    else if (dimension[0] == n2 && dimension[1] == n3 && dimension[2] == n1) {
        swapSizes(0, 1);
        swapSizes(0, 2);
    }
    else if (dimension[0] == n3 && dimension[1] == n1 && dimension[2] == n2) {
        swapSizes(0, 1);
        swapSizes(1, 2);
    }
    else if (dimension[0] == n3 && dimension[1] == n2 && dimension[2] == n1) {
        swapSizes(0, 2);
    }
}

template <int64_t M>
void ModScheme<M>::swapSizes(std::mt19937 &generator) {
    int p1, p2;

    do {
        p1 = ijkDistribution(generator);
        p2 = ijkDistribution(generator);
    } while (p1 == p2);

    swapSizes(p1, p2);
}

template <int64_t M>
void ModScheme<M>::swapSizes(int p1, int p2) {
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
        ModVector u(dimensionNew[0] * dimensionNew[1], M);
        ModVector v(dimensionNew[1] * dimensionNew[2], M);
        ModVector w(dimensionNew[2] * dimensionNew[0], M);

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

    initFlips();
}

template <int64_t M>
void ModScheme<M>::merge(const ModScheme<M> &scheme, int p) {
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
        ModVector u(elementsNew[0], M);
        ModVector v(elementsNew[1], M);
        ModVector w(elementsNew[2], M);

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
        ModVector u(elementsNew[0], M);
        ModVector v(elementsNew[1], M);
        ModVector w(elementsNew[2], M);

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

template <int64_t M>
void ModScheme<M>::project(int p, int q) {
    excludeRow(p, q);
    excludeColumn((p + 2) % 3, q);
    dimension[p]--;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    removeZeroes();
    initFlips();
}

template <int64_t M>
void ModScheme<M>::extend(int p) {
    addRow(p);
    addColumn((p + 2) % 3);

    if (p == 0) {
        for (int i = 0; i < dimension[2]; i++) {
            for (int j = 0; j < dimension[1]; j++) {
                ModVector u((dimension[0] + 1) * dimension[1], dimension[0] * dimension[1] + j, M);
                ModVector v(dimension[1] * dimension[2], j * dimension[2] + i, M);
                ModVector w(dimension[2] * (dimension[0] + 1), i * (dimension[0] + 1) + dimension[0], M);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }
    else if (p == 1) {
        for (int i = 0; i < dimension[0]; i++) {
            for (int j = 0; j < dimension[2]; j++) {
                ModVector u(dimension[0] * (dimension[1] + 1), i * (dimension[1] + 1) + dimension[1], M);
                ModVector v((dimension[1] + 1) * dimension[2], dimension[1] * dimension[2] + j, M);
                ModVector w(dimension[2] * dimension[0], j * dimension[0] + i, M);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }
    else {
        for (int i = 0; i < dimension[0]; i++) {
            for (int j = 0; j < dimension[1]; j++) {
                ModVector u(dimension[0] * dimension[1], i * dimension[1] + j, M);
                ModVector v(dimension[1] * (dimension[2] + 1), j * (dimension[2] + 1) + dimension[2], M);
                ModVector w((dimension[2] + 1) * dimension[0], dimension[2] * dimension[0] + i, M);
                addTriplet(0, 1, 2, u, v, w);
            }
        }
    }

    dimension[p]++;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    initFlips();
}

template <int64_t M>
void ModScheme<M>::product(const ModScheme<M> &scheme2) {
    ModScheme<M> scheme1(*this);

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

                ModVector vector(elements[p], M);

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

template <int64_t M>
void ModScheme<M>::saveJson(const std::string &path) const {
    std::ofstream f(path);

    f << "{" << std::endl;
    f << "    \"n\": [" << dimension[0] << ", " << dimension[1] << ", " << dimension[2] << "]," << std::endl;
    f << "    \"m\": " << rank << "," << std::endl;
    f << "    \"ring\": \"Z" << M << "\"," << std::endl;
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

template <int64_t M>
void ModScheme<M>::saveTxt(const std::string &path) const {
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

template <int64_t M>
void ModScheme<M>::copy(const ModScheme &scheme) {
    rank = scheme.rank;

    for (int i = 0; i < 3; i++) {
        dimension[i] = scheme.dimension[i];
        elements[i] = scheme.elements[i];
        uvw[i].clear();

        for (int index = 0; index < rank; index++) {
            ModVector vector(elements[i], M);

            for (int j = 0; j < elements[i]; j++)
                vector.set(j, scheme.uvw[i][index][j]);

            uvw[i].emplace_back(vector);
        }
    }

    initFlips();
}

template <int64_t M>
bool ModScheme<M>::validate() const {
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                if (!validateEquation(i, j, k))
                    return false;

    return true;
}

template <int64_t M>
bool ModScheme<M>::validateParallel() const {
    bool valid = true;

    #pragma omp parallel for collapse(3) reduction(&&: valid) schedule(dynamic, 32)
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                valid = valid && validateEquation(i, j, k);

    return valid;
}

template <int64_t M>
bool ModScheme<M>::reconstruct(FractionalScheme &scheme) const {
    return false;
}

template <int64_t M>
ModLifter ModScheme<M>::toLift() const {
    std::vector<uint64_t> u(rank * elements[0]);
    std::vector<uint64_t> v(rank * elements[1]);
    std::vector<uint64_t> w(rank * elements[2]);

    for (int index = 0; index < rank; index++) {
        for (int i = 0; i < elements[0]; i++)
            u[index * elements[0] + i] = uvw[0][index][i];

        for (int i = 0; i < elements[1]; i++)
            v[index * elements[1] + i] = uvw[1][index][i];

        for (int i = 0; i < elements[2]; i++)
            w[index * elements[2] + i] = uvw[2][index][i];
    }

    return ModLifter(dimension[0], dimension[1], dimension[2], rank, M, u, v, w, getJakobian());
}

template <int64_t M>
bool ModScheme<M>::lift(FractionalScheme& lifted, int steps) const {
    ModLifter lifter = toLift();
    bool reconstructed = false;

    for (int step = 1; step < steps && !reconstructed && lifter.lift(); step++)
        reconstructed = lifter.reconstruct(lifted) && lifted.validate();

    return reconstructed;
}

template <int64_t M>
bool ModScheme<M>::canLift(int steps) const {
    FractionalScheme lifted;
    return lift(lifted, steps);
}

template <int64_t M>
void ModScheme<M>::initFlips() {
    for (int i = 0; i < 3; i++) {
        flips[i].clear();
        flipsNeg[i].clear();

        for (int index1 = 0; index1 < rank; index1++) {
            for (int index2 = index1 + 1; index2 < rank; index2++) {
                int cmp = uvw[i][index1].compare(uvw[i][index2]);

                if (cmp == 1)
                    flips[i].add(index1, index2);
                else if (cmp == -1)
                    flipsNeg[i].add(index1, index2);
            }
        }
    }
}

template <int64_t M>
void ModScheme<M>::removeZeroes() {
    for (int index = 0; index < rank; index++)
        if (!uvw[0][index] || !uvw[1][index] || !uvw[2][index])
            removeAt(index--);
}

template <int64_t M>
void ModScheme<M>::removeAt(int index) {
    rank--;

    if (index != rank) {
        uvw[0][index] = uvw[0].back();
        uvw[1][index] = uvw[1].back();
        uvw[2][index] = uvw[2].back();
    }

    uvw[0].pop_back();
    uvw[1].pop_back();
    uvw[2].pop_back();
}

template <int64_t M>
void ModScheme<M>::addTriplet(int i, int j, int k, const ModVector &u, const ModVector &v, const ModVector &w) {
    uvw[i].emplace_back(u);
    uvw[j].emplace_back(v);
    uvw[k].emplace_back(w);
    rank++;
}

template <int64_t M>
void ModScheme<M>::excludeColumn(int matrix, int column) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];
    std::vector<int> oldColumns;

    for (int j = 0; j < n2; j++)
        if (j != column)
            oldColumns.push_back(j);

    for (int index = 0; index < rank; index++) {
        ModVector vector(n1 * (n2 - 1), M);

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2 - 1; j++)
                vector.set(i * (n2 - 1) + j, uvw[matrix][index][i * n2 + oldColumns[j]]);

        uvw[matrix][index] = vector;
    }
}

template <int64_t M>
void ModScheme<M>::excludeRow(int matrix, int row) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];
    std::vector<int> oldRows;

    for (int i = 0; i < n1; i++)
        if (i != row)
            oldRows.push_back(i);

    for (int index = 0; index < rank; index++) {
        ModVector vector((n1 - 1) * n2, M);

        for (int i = 0; i < n1 - 1; i++)
            for (int j = 0; j < n2; j++)
                vector.set(i * n2 + j, uvw[matrix][index][oldRows[i] * n2 + j]);

        uvw[matrix][index] = vector;
    }
}

template <int64_t M>
void ModScheme<M>::addColumn(int matrix) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];

    for (int index = 0; index < rank; index++) {
        ModVector vector(n1 * (n2 + 1), M);

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vector.set(i * (n2 + 1) + j, uvw[matrix][index][i * n2 + j]);

        uvw[matrix][index] = vector;
    }
}

template <int64_t M>
void ModScheme<M>::addRow(int matrix) {
    int n1 = dimension[matrix];
    int n2 = dimension[(matrix + 1) % 3];

    for (int index = 0; index < rank; index++) {
        ModVector vector((n1 + 1) * n2, M);

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vector.set(i * n2 + j, uvw[matrix][index][i * n2 + j]);

        uvw[matrix][index] = vector;
    }
}

template <int64_t M>
void ModScheme<M>::flip(int i, int j, int k, int index1, int index2, bool pos) {
    if (pos)
        uvw[j][index1] += uvw[j][index2];
    else
        uvw[j][index1] -= uvw[j][index2];

    uvw[k][index2] -= uvw[k][index1];

    flips[j].remove(index1);
    flips[k].remove(index2);

    flipsNeg[j].remove(index1);
    flipsNeg[k].remove(index2);

    if (!uvw[j][index1] || !uvw[k][index2]) {
        removeZeroes();
        initFlips();
        return;
    }

    for (int index = 0; index < rank; index++) {
        if (index != index1) {
            int cmp = uvw[j][index].compare(uvw[j][index1]);
            if (cmp != 0 && checkFlipReduce(i, k, index, index1, cmp))
                return;

            if (cmp == 1)
                flips[j].add(index1, index);
            else if (cmp == -1)
                flipsNeg[j].add(index1, index);
        }

        if (index != index2) {
            int cmp = uvw[k][index].compare(uvw[k][index2]);
            if (cmp != 0 && checkFlipReduce(i, j, index, index2, cmp))
                return;

            if (cmp == 1)
                flips[k].add(index2, index);
            else if (cmp == -1)
                flipsNeg[k].add(index2, index);
        }
    }
}

template <int64_t M>
void ModScheme<M>::plus(int i, int j, int k, int index1, int index2, int variant) {
    ModVector a1(uvw[i][index1]);
    ModVector b1(uvw[j][index1]);
    ModVector c1(uvw[k][index1]);

    ModVector a2(uvw[i][index2]);
    ModVector b2(uvw[j][index2]);
    ModVector c2(uvw[k][index2]);

    ModVector aAdd = a1 + a2;
    ModVector bAdd = b1 + b2;
    ModVector cAdd = c1 + c2;

    ModVector aSub = a2 - a1;
    ModVector bSub = b2 - b1;
    ModVector cSub = c2 - c1;

    if (variant == 0) {
        uvw[j][index1] = bAdd;
        uvw[i][index2] = aSub;
        addTriplet(i, j, k, a1, b2, cSub);
    }
    else if (variant == 1) {
        uvw[k][index1] = cAdd;
        uvw[j][index2] = bSub;
        addTriplet(i, j, k, aSub, b1, c2);
    }
    else {
        uvw[i][index1] = aAdd;
        uvw[k][index2] = cSub;
        addTriplet(i, j, k, a2, bSub, c1);
    }

    removeZeroes();
    initFlips();
}

template <int64_t M>
void ModScheme<M>::split(int i, int j, int k, int index1, int index2) {
    ModVector u = uvw[i][index1] - uvw[i][index2];
    ModVector v(uvw[j][index1]);
    ModVector w(uvw[k][index1]);

    addTriplet(i, j, k, u, v, w);
    uvw[i][index1] = uvw[i][index2];

    removeZeroes();
    initFlips();
}

template <int64_t M>
bool ModScheme<M>::reduce(int i, int index1, int index2, int sign) {
    if (sign == 0)
        return false;

    if (sign == 1)
        reduceAdd(i, index1, index2);
    else
        reduceSub(i, index1, index2);

    return true;
}

template <int64_t M>
void ModScheme<M>::reduceAdd(int i, int index1, int index2) {
    uvw[i][index1] += uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <int64_t M>
void ModScheme<M>::reduceSub(int i, int index1, int index2) {
    uvw[i][index1] -= uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

template <int64_t M>
bool ModScheme<M>::checkFlipReduce(int i, int j, int index1, int index2, int sign) {
    int cmpI = uvw[i][index1].compare(uvw[i][index2]);
    if (reduce(j, index1, index2, cmpI * sign))
        return true;

    int cmpJ = uvw[j][index1].compare(uvw[j][index2]);
    if (reduce(i, index1, index2, cmpJ * sign))
        return true;

    return false;
}

template <int64_t M>
bool ModScheme<M>::isValidProject(int p, int minN) const {
    return dimension[p] > minN && dimension[(p + 1) % 3] >= minN && dimension[(p + 2) % 3] >= minN;
}

template <int64_t M>
bool ModScheme<M>::isValidExtension(int p, int maxN, int maxRank) const {
    if (rank + dimension[(p + 1) % 3] * dimension[(p + 2) % 3] > maxRank)
        return false;

    int dimensionNew[3] = {dimension[0], dimension[1], dimension[2]};
    dimensionNew[p]++;

    for (int i = 0; i < 3; i++)
        if (dimensionNew[i] > maxN)
            return false;

    return true;
}

template <int64_t M>
bool ModScheme<M>::isValidProduct(const ModScheme<M> &scheme, int maxN, int maxRank) const {
    if (rank * scheme.rank > maxRank)
        return false;

    int dimensionNew[3];

    for (int i = 0; i < 3; i++) {
        dimensionNew[i] = dimension[i] * scheme.dimension[i];

        if (dimensionNew[i] > maxN)
            return false;
    }

    return true;
}

template <int64_t M>
bool ModScheme<M>::isValidMerge(int p, const ModScheme<M> &scheme, int maxN, int maxRank) const {
    int j = (p + 1) % 3;
    int k = (p + 2) % 3;

    int n = dimension[p] + scheme.dimension[p];

    return n <= maxN && dimension[j] == scheme.dimension[j] && dimension[k] == scheme.dimension[k] && rank + scheme.rank <= maxRank;
}

template <int64_t M>
bool ModScheme<M>::validateDimensions() const {
    for (int i = 0; i < 3; i++) {
        if (dimension[i] < 1) {
            std::cout << "Invalid dimension \"" << dimension[i] << "\". Possible dimensions are 1 .." << std::endl;
            return false;
        }

        if (elements[i] < 1) {
            std::cout << "Invalid matrix elements count \"" << elements[i] << "\". Possible counts are 1 .." << std::endl;
            return false;
        }
    }

    if (rank < 1) {
        std::cout << "Invalid rank \"" << rank << "\"" << std::endl;
        return false;
    }

    return true;
}

template <int64_t M>
bool ModScheme<M>::validateEquation(int i, int j, int k) const {
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

    return equation % M == target;
}

template <int64_t M>
void ModScheme<M>::saveMatrix(std::ofstream &f, std::string name, const std::vector<ModVector> &vectors) const {
    f << "    \"" << name << "\": [" << std::endl;

    for (size_t index = 0; index < vectors.size(); index++)
        f << "        [" << vectors[index] << "]" << (index < vectors.size() - 1 ? "," : "") << std::endl;

    f << "    ]";
}

template <int64_t M>
ModMatrix ModScheme<M>::getJakobian() const {
    int rows = elements[0] * elements[1] * elements[2];
    int columns = rank * (elements[0] + elements[1] + elements[2]);
    ModMatrix jakobian(rows, columns, M);

    int vOffset = elements[0] * rank;
    int wOffset = (elements[0] + elements[1]) * rank;

    for (int i = 0; i < elements[0]; i++) {
        for (int j = 0; j < elements[1]; j++) {
            for (int k = 0; k < elements[2]; k++) {
                int row = (i * elements[1] + j) * elements[2] + k;

                for (int index = 0; index < rank; index++) {
                    int u = uvw[0][index][i];
                    int v = uvw[1][index][j];
                    int w = uvw[2][index][k];

                    jakobian(row, i * rank + index) = v * w;
                    jakobian(row, vOffset + j * rank + index) = u * w;
                    jakobian(row, wOffset + k * rank + index) = u * v;
                }
            }
        }
    }

    return jakobian;
}
