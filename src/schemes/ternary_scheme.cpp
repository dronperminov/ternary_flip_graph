#include "ternary_scheme.h"

TernaryScheme::TernaryScheme() : boolDistribution(0, 1), ijkDistribution(0, 2) {

}

TernaryScheme::TernaryScheme(const TernaryScheme &scheme) {
    copy(scheme);
}

void TernaryScheme::initializeNaive(int n1, int n2, int n3) {
    dimension[0] = n1;
    dimension[1] = n2;
    dimension[2] = n3;

    elements[0] = n1 * n2;
    elements[1] = n2 * n3;
    elements[2] = n3 * n1;

    rank = n1 * n2 * n3;

    uvw[0].clear();
    uvw[1].clear();
    uvw[2].clear();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int k = 0; k < n2; k++) {
                uvw[0].emplace_back(TernaryVector<vec_type>(n1 * n2, i * n2 + k));
                uvw[1].emplace_back(TernaryVector<vec_type>(n2 * n3, k * n3 + j));
                uvw[2].emplace_back(TernaryVector<vec_type>(n3 * n1, j * n1 + i));
            }
        }
    }

    initFlips();
}

bool TernaryScheme::read(const std::string &path) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    f >> dimension[0] >> dimension[1] >> dimension[2] >> rank;

    int maxSize = sizeof(vec_type) * 8;

    for (int i = 0; i < 3; i++) {
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

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

    for (int i = 0; i < 3; i++) {
        for (int index = 0; index < rank; index++) {
            TernaryVector<vec_type> vector(elements[i]);
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

int TernaryScheme::getRank() const {
    return rank;
}

int TernaryScheme::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += uvw[i][index].nonZeroCount();

    return count - 2 * rank - elements[2];
}

int TernaryScheme::getDimension(int index) const {
    return dimension[index];
}

int TernaryScheme::getAvailableFlips() const {
    return flips[0].size() + flips[1].size() + flips[2].size();
}

bool TernaryScheme::tryFlip(std::mt19937 &generator) {
    size_t size = flips[0].size() + flips[1].size() + flips[2].size();

    if (!size)
        return false;

    std::vector<int> indices(size);

    for (size_t i = 0; i < size; i++)
        indices[i] = i;

    std::shuffle(indices.begin(), indices.end(), generator);

    for (size_t p = 0; p < size; p++) {
        size_t index = indices[p];
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

bool TernaryScheme::tryPlus(std::mt19937 &generator) {
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

bool TernaryScheme::tryReduce() {
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

bool TernaryScheme::validate() const {
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

void TernaryScheme::saveJson(const std::string &path) const {
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

void TernaryScheme::saveTxt(const std::string &path) const {
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

void TernaryScheme::copy(const TernaryScheme &scheme) {
    rank = scheme.rank;

    for (int i = 0; i < 3; i++) {
        dimension[i] = scheme.dimension[i];
        elements[i] = scheme.elements[i];
        uvw[i].clear();

        for (int index = 0; index < rank; index++) {
            TernaryVector<vec_type> vector(elements[i]);

            for (int j = 0; j < elements[i]; j++)
                vector.set(j, scheme.uvw[i][index][j]);

            uvw[i].emplace_back(vector);
        }
    }

    initFlips();
}

void TernaryScheme::initFlips() {
    for (int i = 0; i < 3; i++) {
        flips[i].clear();

        for (int index1 = 0; index1 < rank; index1++)
            for (int index2 = index1 + 1; index2 < rank; index2++)
                if (uvw[i][index1] == uvw[i][index2])
                    flips[i].add(index1, index2);
    }
}

void TernaryScheme::removeZeroes() {
    for (int index = 0; index < rank; index++)
        if (!uvw[0][index] || !uvw[1][index] || !uvw[2][index])
            removeAt(index--);
}

void TernaryScheme::removeAt(int index) {
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

void TernaryScheme::addTriplet(int i, int j, int k, const TernaryVector<vec_type> &u, const TernaryVector<vec_type> &v, const TernaryVector<vec_type> &w) {
    uvw[i].emplace_back(u);
    uvw[j].emplace_back(v);
    uvw[k].emplace_back(w);
    rank++;
}

void TernaryScheme::flip(int i, int j, int k, int index1, int index2) {
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

bool TernaryScheme::plus(int i, int j, int k, int index1, int index2, int variant) {
    TernaryVector<vec_type> a1(uvw[i][index1]);
    TernaryVector<vec_type> b1(uvw[j][index1]);
    TernaryVector<vec_type> c1(uvw[k][index1]);

    TernaryVector<vec_type> a2(uvw[i][index2]);
    TernaryVector<vec_type> b2(uvw[j][index2]);
    TernaryVector<vec_type> c2(uvw[k][index2]);

    TernaryVector<vec_type> aAdd = a1 + a2;
    TernaryVector<vec_type> bAdd = b1 + b2;
    TernaryVector<vec_type> cAdd = c1 + c2;

    TernaryVector<vec_type> aSub = a2 - a1;
    TernaryVector<vec_type> bSub = b2 - b1;
    TernaryVector<vec_type> cSub = c2 - c1;

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

void TernaryScheme::reduceAdd(int i, int index1, int index2) {
    uvw[i][index1] += uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

void TernaryScheme::reduceSub(int i, int index1, int index2) {
    uvw[i][index1] -= uvw[i][index2];
    bool isZero = !uvw[i][index1];

    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

bool TernaryScheme::checkFlipReduce(int i, int j, int index1, int index2) {
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

bool TernaryScheme::fixSigns() {
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

bool TernaryScheme::validateEquation(int i, int j, int k) const {
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

void TernaryScheme::saveMatrix(std::ofstream &f, std::string name, const std::vector<TernaryVector<vec_type>> &vectors) const {
    f << "    \"" << name << "\": [" << std::endl;

    for (size_t index = 0; index < vectors.size(); index++)
        f << "        [" << vectors[index] << "]" << (index < vectors.size() - 1 ? "," : "") << std::endl;

    f << "    ]";
}
