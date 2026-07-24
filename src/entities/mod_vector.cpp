#include "mod_vector.h"

ModVector::ModVector(int64_t mod) {
    this->n = 0;
    this->mod = mod;
}

ModVector::ModVector(int n, int64_t mod) {
    this->n = n;
    this->mod = mod;
    this->values = std::vector<int64_t>(n, 0);
}

ModVector::ModVector(int n, int index, int64_t mod) {
    this->n = n;
    this->mod = mod;
    this->values = std::vector<int64_t>(n, 0);
    this->values[index] = 1;
}

ModVector::ModVector(int n, int *values, int64_t mod) {
    this->n = n;
    this->mod = mod;
    this->values = std::vector<int64_t>(n, 0);

    for (int i = 0; i < n; i++)
        set(i, values[i]);
}

ModVector::ModVector(const ModVector& vector) {
    n = vector.n;
    mod = vector.mod;
    values = vector.values;
}

void ModVector::set(int index, int value) {
    value %= mod;
    if (value < 0)
        value += mod;

    values[index] = value;
}

int ModVector::compare(const ModVector &vector) const {
    bool eq = true;
    bool neg = true;

    for (int i = 0; i < n && (eq || neg); i++) {
        if (values[i] != vector.values[i])
            eq = false;

        if (values[i] != (mod - vector.values[i]) % mod)
            neg = false;
    }

    if (eq)
        return 1;

    if (neg)
        return -1;

    return 0;
}

bool ModVector::operator==(const ModVector &vector) const {
    for (int i = 0; i < n; i++)
        if (values[i] != vector.values[i])
            return false;

    return true;
}

bool ModVector::operator!=(const ModVector &vector) const {
    for (int i = 0; i < n; i++)
        if (values[i] != vector.values[i])
            return true;

    return false;
}

int ModVector::operator[](int index) const {
    return values[index];
}

ModVector ModVector::operator+(const ModVector &vector) const {
    ModVector result(n, mod);

    for (int i = 0; i < n; i++)
        result.values[i] = (values[i] + vector.values[i]) % mod;

    return result;
}

ModVector ModVector::operator-(const ModVector &vector) const {
    ModVector result(n, mod);

    for (int i = 0; i < n; i++)
        result.values[i] = (values[i] + mod - vector.values[i]) % mod;

    return result;
}

ModVector ModVector::operator*(int value) const {
    ModVector result(n, mod);
    value %= mod;

    for (int i = 0; i < n; i++)
        result.values[i] = (values[i] * value) % mod;

    return result;
}

ModVector ModVector::operator-() const {
    ModVector result(n, mod);

    for (int i = 0; i < n; i++)
        result.values[i] = (mod - values[i]) % mod;

    return result;
}

ModVector& ModVector::operator+=(const ModVector &vector) {
    for (int i = 0; i < n; i++)
        values[i] = (values[i] + vector.values[i]) % mod;

    return *this;
}

ModVector& ModVector::operator-=(const ModVector &vector) {
    for (int i = 0; i < n; i++)
        values[i] = (values[i] + mod - vector.values[i]) % mod;

    return *this;
}

ModVector& ModVector::operator*=(int value) {
    value %= mod;

    for (int i = 0; i < n; i++)
        values[i] = (values[i] * value) % mod;

    return *this;
}

int ModVector::nonZeroCount() const {
    int count = 0;

    for (int i = 0; i < n; i++)
        if (values[i])
            count++;

    return count;
}

ModVector::operator bool() const {
    for (int i = 0; i < n; i++)
        if (values[i])
            return true;

    return false;
}

std::ostream& operator<<(std::ostream &os, const ModVector &vector) {
    for (int i = 0; i < vector.n; i++) {
        if (i > 0)
            os << ", ";

        os << vector.values[i];
    }

    return os;
}

std::istream& operator>>(std::istream &is, ModVector &vector) {
    int value;

    for (int i = 0; i < vector.n; i++) {
        is >> value;
        vector.set(i, value);
    }

    return is;
}
