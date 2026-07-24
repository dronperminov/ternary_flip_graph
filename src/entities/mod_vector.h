#pragma once

#include <iostream>
#include <vector>
#include <cstdint>

struct ModVector {
    int mod;
    int n;
    std::vector<int64_t> values;

    ModVector(int64_t mod);
    ModVector(int n, int64_t mod);
    ModVector(int n, int index, int64_t mod);
    ModVector(int n, int *values, int64_t mod);
    ModVector(const ModVector& vector);

    void set(int index, int value);

    int compare(const ModVector &vector) const;
    bool operator==(const ModVector &vector) const;
    bool operator!=(const ModVector &vector) const;
    int operator[](int index) const;

    ModVector operator+(const ModVector &vector) const;
    ModVector operator-(const ModVector &vector) const;
    ModVector operator*(int value) const;
    ModVector operator-() const;

    ModVector& operator+=(const ModVector &vector);
    ModVector& operator-=(const ModVector &vector);
    ModVector& operator*=(int value);

    int nonZeroCount() const;

    explicit operator bool() const;

    friend std::ostream& operator<<(std::ostream &os, const ModVector &vector);
    friend std::istream& operator>>(std::istream &is, ModVector &vector);
};
