#include "fractional_scheme.h"

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

bool FractionalScheme::validateParallel() const {
    bool valid = true;

    #pragma omp parallel for collapse(3) reduction(&&: valid) schedule(dynamic, 32)
    for (int i = 0; i < elements[0]; i++)
        for (int j = 0; j < elements[1]; j++)
            for (int k = 0; k < elements[2]; k++)
                valid = valid && validateEquation(i, j, k);

    return valid;
}

bool FractionalScheme::read(const std::string &path, bool checkCorrectness, bool integer) {
    std::ifstream f(path);

    if (!f) {
        std::cout << "Unable open file \"" << path << "\"" << std::endl;
        return false;
    }

    bool valid = read(f, checkCorrectness, integer);
    f.close();

    if (!valid) {
        std::cout << "Invalid scheme in the file \"" << path << "\"" << std::endl;
        return false;
    }

    return true;
}

bool FractionalScheme::read(std::istream &is, bool checkCorrectness, bool integer) {
    is >> dimension[0] >> dimension[1] >> dimension[2] >> rank;

    for (int i = 0; i < 3; i++)
        elements[i] = dimension[i] * dimension[(i + 1) % 3];

    int numerator;
    int denominator = 1;

    for (int i = 0; i < 3; i++) {
        uvw[i].resize(rank * elements[i]);

        for (int j = 0; j < rank * elements[i]; j++) {
            is >> numerator;
            if (!integer)
                is >> denominator;

            uvw[i][j] = Fraction(numerator, denominator);
        }
    }

    if (checkCorrectness && !validateParallel())
        return false;

    initFlips();
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

int FractionalScheme::getAvailableFlips() const {
    size_t sizePos = flips[0].size() + flips[1].size() + flips[2].size();
    size_t sizeNeg = flipsNeg[0].size() + flipsNeg[1].size() + flipsNeg[2].size();
    return sizePos + sizeNeg;
}

int FractionalScheme::getAvailableFlips(int index) const {
    return flips[index].size() + flipsNeg[index].size();
}

int FractionalScheme::getFractionsCount() const {
    return getFractionsCount(0) + getFractionsCount(1) + getFractionsCount(2);
}

int FractionalScheme::getFractionsCount(int index) const {
    int count = 0;

    for (int i = 0; i < rank * elements[index]; i++)
        if (!uvw[index][i].isInteger())
            count++;

    return count;
}

int FractionalScheme::getComplexity() const {
    int complexity = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            if (uvw[i][j] != 0)
                complexity++;

    return complexity - 2 * rank - elements[2];
}

int64_t FractionalScheme::getWeight() const {
    int64_t weight = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            weight += abs(uvw[i][j].numerator()) * uvw[i][j].denominator();

    return weight;
}

int FractionalScheme::getMaxAbsNumerator() const {
    int maxAbsNumerator = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            maxAbsNumerator = std::max(maxAbsNumerator, abs(uvw[i][j].numerator()));

    return maxAbsNumerator;
}

int FractionalScheme::getAbsNumeratorCount(int numerator) const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            if (abs(uvw[i][j].numerator()) == numerator)
                count++;

    return count;
}

int FractionalScheme::getMaxDenominator() const {
    int denominator = 1;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            denominator = std::max(denominator, uvw[i][j].denominator());

    return denominator;
}

std::string FractionalScheme::getRing() const {
    if (isTernary())
        return "ZT";

    if (isInteger())
        return "Z";

    return "Q";
}

std::string FractionalScheme::getHash() const {
    std::vector<std::string> lines;

    for (int index = 0; index < rank; index++) {
        std::stringstream ss;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < elements[i]; j++)
                ss << uvw[i][index * elements[i] + j].pretty();

        lines.push_back(ss.str());
    }

    std::sort(lines.begin(), lines.end());
    std::stringstream hash;
    for (int index = 0; index < rank; index++)
        hash << lines[index];

    return hash.str();
}

std::string FractionalScheme::getUniqueValues() const {
    std::unordered_set<Fraction> unique;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < rank * elements[i]; j++)
            unique.insert(uvw[i][j]);

    std::vector<Fraction> uniqueValues(unique.begin(), unique.end());
    std::sort(uniqueValues.begin(), uniqueValues.end());
    std::stringstream ss;

    ss << "{";
    for (size_t i = 0; i < uniqueValues.size(); i++)
        ss << (i > 0 ? ", " : "") << uniqueValues[i].pretty();
    ss << "}";

    return ss.str();
}

double FractionalScheme::getFrobeniusNorm() const {
    double norm = 0;

    for (int index = 0; index < rank; index++) {
        Matrix u(dimension[0], dimension[1]);
        Matrix v(dimension[1], dimension[2]);
        Matrix w(dimension[2], dimension[0]);

        for (int i = 0; i < elements[0]; i++)
            u[i] = uvw[0][index * elements[0] + i];

        for (int i = 0; i < elements[1]; i++)
            v[i] = uvw[1][index * elements[1] + i];

        for (int i = 0; i < elements[2]; i++)
            w[i] = uvw[2][index * elements[2] + i];

        Matrix uut = u * u.transpose();
        Matrix vvt = v * v.transpose();
        Matrix wwt = w * w.transpose();

        norm += std::sqrt((uut.trace() * vvt.trace() * wwt.trace()).toDouble());
    }

    return norm;
}

FlipStructure FractionalScheme::getOptimalStructure(std::mt19937 &generator, int iterations, double eps) const {
    FlipStructureOptimizer optimizer(dimension[0], dimension[1], dimension[2], rank);

    for (int i = 0; i < 3; i++) {
        for (size_t j = 0; j < flips[i].size(); j++)
            optimizer.add(i, flips[i].index1(j), flips[i].index2(j));

        for (size_t j = 0; j < flipsNeg[i].size(); j++)
            optimizer.add(i, flipsNeg[i].index1(j), flipsNeg[i].index2(j));
    }

    return optimizer.optimize(generator, iterations, eps);
}

std::string FractionalScheme::getTypeInvariant() const {
    std::vector<Ranks> ranks;

    for (int index = 0; index < rank; index++) {
        Matrix u(dimension[0], dimension[1]);
        Matrix v(dimension[1], dimension[2]);
        Matrix w(dimension[2], dimension[0]);

        for (int i = 0; i < elements[0]; i++)
            u[i] = uvw[0][index * elements[0] + i];

        for (int i = 0; i < elements[1]; i++)
            v[i] = uvw[1][index * elements[1] + i];

        for (int i = 0; i < elements[2]; i++)
            w[i] = uvw[2][index * elements[2] + i];

        ranks.push_back({u.rank(), v.rank(), w.rank()});
    }

    std::unordered_map<Ranks, int> powers;
    for (int index = 0; index < rank; index++) {
        auto result = powers.find(ranks[index]);

        if (result == powers.end())
            powers[ranks[index]] = 1;
        else
            result->second++;
    }

    std::vector<std::string> types;
    for (const auto &power: powers) {
        std::stringstream type;

        if (power.second > 1)
            type << power.second;

        type << power.first;
        types.push_back(type.str());
    }

    std::sort(types.begin(), types.end());
    std::stringstream type;

    for (size_t i = 0; i < types.size(); i++) {
        if (i > 0)
            type << "+";

        type << types[i];
    }

    return type.str();
}

bool FractionalScheme::tryFlip(std::mt19937 &generator) {
    size_t sizePos = flips[0].size() + flips[1].size() + flips[2].size();
    size_t sizeNeg = flipsNeg[0].size() + flipsNeg[1].size() + flipsNeg[2].size();
    size_t size = sizePos + sizeNeg;

    if (!size)
        return false;

    size_t index = generator() % size;
    int i, j, k, index1, index2;

    if (index < sizePos) {
        selectFlip(flips, index, i, j, k, index1, index2, generator);
        flip(i, j, k, index1, index2, false);
    }
    else {
        selectFlip(flipsNeg, index - sizePos, i, j, k, index1, index2, generator);
        flip(i, j, k, index1, index2, true);
    }

    return true;
}

void FractionalScheme::plus(std::mt19937 &generator) {
    int index1 = generator() % rank;
    int index2 = generator() % rank;

    while (index1 == index2 || isEqualMatrices(0, index1, index2) || isEqualMatrices(1, index1, index2) || isEqualMatrices(2, index1, index2)) {
        index1 = generator() % rank;
        index2 = generator() % rank;
    }

    int permutation[3] = {0, 1, 2};
    std::shuffle(permutation, permutation + 3, generator);

    plus(permutation[0], permutation[1], permutation[2], index1, index2, generator() % 3);
}

void FractionalScheme::split(std::mt19937 &generator, const std::vector<Fraction> &values) {
    int permutation[3] = {0, 1, 2};
    std::shuffle(permutation, permutation + 3, generator);

    int i = permutation[0];
    int j = permutation[1];
    int k = permutation[2];
    int index = generator() % rank;

    std::vector<Fraction> u1(elements[i]), u2(elements[i]);
    std::vector<Fraction> v(uvw[j].begin() + index * elements[j], uvw[j].begin() + (index + 1) * elements[j]);
    std::vector<Fraction> w(uvw[k].begin() + index * elements[k], uvw[k].begin() + (index + 1) * elements[k]);

    for (int ind = 0; ind < elements[i]; ind++) {
        Fraction value = values[generator() % values.size()];

        if (boolDistribution(generator)) {
            u1[ind] = value;
            u2[ind] = uvw[i][index * elements[i] + ind] - value;
        }
        else {
            u1[ind] = uvw[i][index * elements[i] + ind] - value;
            u2[ind] = value;
        }
    }

    std::copy(u1.begin(), u1.end(), uvw[i].begin() + index * elements[i]);
    addTriplet(i, j, k, u2, v, w);

    removeZeroes();
    initFlips();
}

void FractionalScheme::sandwiching(const Matrix &u, const Matrix &v, const Matrix &w, const Matrix &u1, const Matrix &v1, const Matrix &w1) {
    for (int index = 0; index < rank; index++) {
        Matrix ui(dimension[0], dimension[1]);
        Matrix vi(dimension[1], dimension[2]);
        Matrix wi(dimension[2], dimension[0]);

        for (int i = 0; i < elements[0]; i++)
            ui[i] = uvw[0][index * elements[0] + i];

        for (int i = 0; i < elements[1]; i++)
            vi[i] = uvw[1][index * elements[1] + i];

        for (int i = 0; i < elements[2]; i++)
            wi[i] = uvw[2][index * elements[2] + i];

        ui = u * ui * v1;
        vi = v * vi * w1;
        wi = w * wi * u1;

        for (int i = 0; i < elements[0]; i++)
            uvw[0][index * elements[0] + i] = ui[i];

        for (int i = 0; i < elements[1]; i++)
            uvw[1][index * elements[1] + i] = vi[i];

        for (int i = 0; i < elements[2]; i++)
            uvw[2][index * elements[2] + i] = wi[i];
    }
}

void FractionalScheme::scale(int index, const Fraction &alpha, const Fraction &beta, const Fraction &gamma) {
    Fraction scale[3] = {alpha, beta, gamma};

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < elements[i]; j++)
            uvw[i][index * elements[i] + j] *= scale[i];
}

void FractionalScheme::fixFractions() {
    for (int index = 0; index < rank; index++) {
        int64_t lcm[3] = {1, 1, 1};
        int64_t gcd[3] = {0, 0, 0};

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < elements[i]; j++) {
                lcm[i] = std::lcm(lcm[i], uvw[i][index * elements[i] + j].denominator());

                int num = uvw[i][index * elements[i] + j].numerator();
                if (num)
                    gcd[i] = std::gcd(gcd[i], std::abs(num));
            }
        }

        int64_t gcdP = gcd[0] * gcd[1] * gcd[2];
        int64_t lcmP = lcm[0] * lcm[1] * lcm[2];
        if (lcmP == 1 || gcdP % lcmP != 0)
            continue;

        Fraction alpha(lcm[0], std::gcd(lcmP, gcd[0]));
        Fraction beta(lcm[1], std::gcd(lcmP, gcd[1]));
        Fraction gamma(lcm[2], std::gcd(lcmP, gcd[2]));

        scale(index, alpha, beta, gamma);
    }
}

void FractionalScheme::copy(const FractionalScheme &scheme) {
    rank = scheme.rank;

    for (int i = 0; i < 3; i++) {
        dimension[i] = scheme.dimension[i];
        elements[i] = scheme.elements[i];
        uvw[i].assign(rank * elements[i], 0);

        for (int j = 0; j < rank * elements[i]; j++)
            uvw[i][j] = scheme.uvw[i][j];
    }

    initFlips();
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

void FractionalScheme::saveJson(const std::string &path, bool withInvariants) const {
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

    if (withInvariants) {
        f << "," << std::endl << "    \"type\": \"" << getTypeInvariant() << "\"";
    }

    f << std::endl;
    f << "}" << std::endl;

    f.close();
}

void FractionalScheme::saveTxt(const std::string &path) const {
    std::ofstream f(path);
    f << dimension[0] << " " << dimension[1] << " " << dimension[2] << " " << rank << std::endl;

    bool fractional = !isInteger();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < rank * elements[i]; j++) {
            if (j > 0)
                f << " ";

            f << uvw[i][j].numerator();

            if (fractional)
                f << " " << uvw[i][j].denominator();
        }

        f << std::endl;
    }

    f.close();
}

void FractionalScheme::save(const std::string &path) const {
    if (path.length() >= 5 && !path.compare(path.length() - 5, 5, ".json")) {
        saveJson(path);
        return;
    }

    saveTxt(path);
}

void FractionalScheme::initFlips() {
    for (int i = 0; i < 3; i++) {
        flips[i].clear();
        flipsNeg[i].clear();

        for (int index1 = 0; index1 < rank; index1++) {
            for (int index2 = index1 + 1; index2 < rank; index2++) {
                if (isEqualMatrices(i, index1, index2)) {
                    flips[i].add(index1, index2);
                }
                else if (isInverseMatrices(i, index1, index2)) {
                    flipsNeg[i].add(index1, index2);
                }
            }
        }
    }
}

void FractionalScheme::removeZeroes() {
    for (int index = 0; index < rank; index++)
        if (isZeroMatrix(0, index) || isZeroMatrix(1, index) || isZeroMatrix(2, index))
            removeAt(index--);
}

void FractionalScheme::removeAt(int index) {
    rank--;

    if (index != rank) {
        for (int i = 0; i < 3; i++) {
            int curr = index * elements[i];
            int last = rank * elements[i];

            for (int j = 0; j < elements[i]; j++)
                uvw[i][curr + j] = uvw[i][last + j];
        }
    }

    for (int i = 0; i < 3; i++)
        uvw[i].resize(rank * elements[i]);
}

void FractionalScheme::addTriplet(int i, int j, int k, const std::vector<Fraction> &u, const std::vector<Fraction> &v, const std::vector<Fraction> &w) {
    uvw[i].insert(uvw[i].end(), u.begin(), u.end());
    uvw[j].insert(uvw[j].end(), v.begin(), v.end());
    uvw[k].insert(uvw[k].end(), w.begin(), w.end());
    rank++;
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

bool FractionalScheme::isEqualMatrices(int p, int index1, int index2) const {
    for (int i = 0; i < elements[p]; i++)
        if (uvw[p][index1 * elements[p] + i] != uvw[p][index2 * elements[p] + i])
            return false;

    return true;
}

bool FractionalScheme::isInverseMatrices(int p, int index1, int index2) const {
    for (int i = 0; i < elements[p]; i++)
        if (uvw[p][index1 * elements[p] + i] != -uvw[p][index2 * elements[p] + i])
            return false;

    return true;
}

int FractionalScheme::compareMatrices(int p, int index1, int index2) const {
    if (isEqualMatrices(p, index1, index2))
        return 1;

    if (isInverseMatrices(p, index1, index2))
        return -1;

    return 0;
}

bool FractionalScheme::isZeroMatrix(int p, int index) const {
    for (int i = 0; i < elements[p]; i++)
        if (uvw[p][index * elements[p] + i])
            return false;

    return true;
}

void FractionalScheme::selectFlip(FlipSet *flips, size_t index, int &i, int &j, int &k, int &index1, int &index2, std::mt19937 &generator) {
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

    index1 = flips[i].index1(index);
    index2 = flips[i].index2(index);

    if (boolDistribution(generator))
        std::swap(j, k);

    if (boolDistribution(generator))
        std::swap(index1, index2);
}

void FractionalScheme::flip(int i, int j, int k, int index1, int index2, bool inverse) {
    if (inverse) {
        for (int index = 0; index < elements[j]; index++)
            uvw[j][index1 * elements[j] + index] += uvw[j][index2 * elements[j] + index];
    }
    else {
        for (int index = 0; index < elements[j]; index++)
            uvw[j][index1 * elements[j] + index] -= uvw[j][index2 * elements[j] + index];
    }

    for (int index = 0; index < elements[k]; index++)
        uvw[k][index2 * elements[k] + index] += uvw[k][index1 * elements[k] + index];

    flips[j].remove(index1);
    flips[k].remove(index2);
    flipsNeg[j].remove(index1);
    flipsNeg[k].remove(index2);

    if (isZeroMatrix(j, index1) || isZeroMatrix(k, index2)) {
        removeZeroes();
        initFlips();
        return;
    }

    for (int index = 0; index < rank; index++) {
        if (index != index1) {
            int cmp = compareMatrices(j, index, index1);
            if (cmp != 0 && checkFlipReduce(i, k, index, index1, cmp))
                return;

            if (cmp == 1) {
                flips[j].add(index1, index);
            }
            else if (cmp == -1) {
                flipsNeg[j].add(index1, index);
            }
        }

        if (index != index2) {
            int cmp = compareMatrices(k, index, index2);
            if (cmp != 0 && checkFlipReduce(i, j, index, index2, cmp))
                return;

            if (cmp == 1) {
                flips[k].add(index2, index);
            }
            else if (cmp == -1) {
                flipsNeg[k].add(index2, index);
            }
        }
    }
}

void FractionalScheme::plus(int i, int j, int k, int index1, int index2, int variant) {
    std::vector<Fraction> a1(uvw[i].begin() + index1 * elements[i], uvw[i].begin() + (index1 + 1) * elements[i]);
    std::vector<Fraction> b1(uvw[j].begin() + index1 * elements[j], uvw[j].begin() + (index1 + 1) * elements[j]);
    std::vector<Fraction> c1(uvw[k].begin() + index1 * elements[k], uvw[k].begin() + (index1 + 1) * elements[k]);

    std::vector<Fraction> a2(uvw[i].begin() + index2 * elements[i], uvw[i].begin() + (index2 + 1) * elements[i]);
    std::vector<Fraction> b2(uvw[j].begin() + index2 * elements[j], uvw[j].begin() + (index2 + 1) * elements[j]);
    std::vector<Fraction> c2(uvw[k].begin() + index2 * elements[k], uvw[k].begin() + (index2 + 1) * elements[k]);

    std::vector<Fraction> aAdd = addVectors(a1, a2);
    std::vector<Fraction> bAdd = addVectors(b1, b2);
    std::vector<Fraction> cAdd = addVectors(c1, c2);

    std::vector<Fraction> aSub = subVectors(a2, a1);
    std::vector<Fraction> bSub = subVectors(b2, b1);
    std::vector<Fraction> cSub = subVectors(c2, c1);

    if (variant == 0) {
        std::copy(bAdd.begin(), bAdd.end(), uvw[j].begin() + index1 * elements[j]);
        std::copy(aSub.begin(), aSub.end(), uvw[i].begin() + index2 * elements[i]);
        addTriplet(i, j, k, a1, b2, cSub);
    }
    else if (variant == 1) {
        std::copy(cAdd.begin(), cAdd.end(), uvw[k].begin() + index1 * elements[k]);
        std::copy(bSub.begin(), bSub.end(), uvw[j].begin() + index2 * elements[j]);
        addTriplet(i, j, k, aSub, b1, c2);
    }
    else {
        std::copy(aAdd.begin(), aAdd.end(), uvw[i].begin() + index1 * elements[i]);
        std::copy(cSub.begin(), cSub.end(), uvw[k].begin() + index2 * elements[k]);
        addTriplet(i, j, k, a2, bSub, c1);
    }

    removeZeroes();
    initFlips();
}

bool FractionalScheme::reduce(int i, int index1, int index2, int sign) {
    if (sign == 0)
        return false;

    if (sign == 1)
        reduceAdd(i, index1, index2);
    else
        reduceSub(i, index1, index2);

    return true;
}

void FractionalScheme::reduceAdd(int i, int index1, int index2) {
    for (int index = 0; index < elements[i]; index++)
        uvw[i][index1 * elements[i] + index] += uvw[i][index2 * elements[i] + index];

    bool isZero = isZeroMatrix(i, index1);
    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

void FractionalScheme::reduceSub(int i, int index1, int index2) {
    for (int index = 0; index < elements[i]; index++)
        uvw[i][index1 * elements[i] + index] -= uvw[i][index2 * elements[i] + index];

    bool isZero = isZeroMatrix(i, index1);
    removeAt(index2);

    if (isZero)
        removeZeroes();

    initFlips();
}

bool FractionalScheme::checkFlipReduce(int i, int j, int index1, int index2, int sign) {
    int cmpI = compareMatrices(i, index1, index2);
    if (reduce(j, index1, index2, cmpI * sign))
        return true;

    int cmpJ = compareMatrices(j, index1, index2);
    if (reduce(i, index1, index2, cmpJ * sign))
        return true;

    return false;
}

std::vector<Fraction> FractionalScheme::addVectors(const std::vector<Fraction> &a, const std::vector<Fraction> &b) const {
    std::vector<Fraction> result(a.size());

    for (size_t i = 0; i < a.size(); i++)
        result[i] = a[i] + b[i];

    return result;
}

std::vector<Fraction> FractionalScheme::subVectors(const std::vector<Fraction> &a, const std::vector<Fraction> &b) const {
    std::vector<Fraction> result(a.size());

    for (size_t i = 0; i < a.size(); i++)
        result[i] = a[i] - b[i];

    return result;
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
