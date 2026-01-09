#pragma once

#include <iostream>
#include <sstream>
#include <numeric>
#include <string>
#include <stdexcept>
#include <cassert>

class Fraction {
    int64_t num;
    int64_t den;
public:
    Fraction(int64_t numerator = 0, int64_t denominator = 1);

    int numerator() const;
    int denominator() const;
    std::string pretty() const;

    bool isInteger() const;
    bool isTernaryInteger() const;
    bool reconstruct(int64_t a, int64_t mod, int64_t bound);

    Fraction operator-() const;
    Fraction operator+(const Fraction &fraction) const;
    Fraction operator-(const Fraction &fraction) const;
    Fraction operator*(const Fraction &fraction) const;
    Fraction operator/(const Fraction &fraction) const;

    Fraction& operator+=(const Fraction &fraction);
    Fraction& operator-=(const Fraction &fraction);
    Fraction& operator*=(const Fraction &fraction);
    Fraction& operator/=(const Fraction &fraction);

    bool operator==(const Fraction &fraction) const;
    bool operator==(int value) const;

    bool operator!=(const Fraction &fraction) const;
    bool operator!=(int value) const;

    bool operator>(const Fraction& fraction) const;
    bool operator<(const Fraction& fraction) const;

    friend std::ostream& operator<<(std::ostream &os, const Fraction &fraction);
    friend std::istream& operator>>(std::istream &is, Fraction &fraction);
private:
    void normalize();
};

Fraction abs(const Fraction &fraction);
