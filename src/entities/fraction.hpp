#pragma once

#include <iostream>
#include <numeric>
#include <stdexcept>

class Fraction {
    int64_t num;
    int64_t den;
public:
    Fraction(int64_t numerator = 0, int64_t denominator = 1);

    int numerator() const;
    int denominator() const;

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

    friend std::ostream& operator<<(std::ostream &os, const Fraction &fraction);
private:
    void normalize();
};

Fraction::Fraction(int64_t numerator, int64_t denominator) {
    num = denominator > 0 ? numerator : -numerator;
    den = abs(denominator);
    normalize();
}

int Fraction::numerator() const {
    return num;
}

int Fraction::denominator() const {
    return den;
}

Fraction Fraction::operator-() const {
    return Fraction(-num, den);
}

Fraction Fraction::operator+(const Fraction &fraction) const {
    int64_t numerator = num * fraction.den + fraction.num * den;
    int64_t denominator = den * fraction.den;
    return Fraction(numerator, denominator);
}

Fraction Fraction::operator-(const Fraction &fraction) const {
    int64_t numerator = num * fraction.den - fraction.num * den;
    int64_t denominator = den * fraction.den;
    return Fraction(numerator, denominator);
}

Fraction Fraction::operator*(const Fraction &fraction) const {
    int64_t gcd1 = std::gcd(abs(num), fraction.den);
    int64_t gcd2 = std::gcd(abs(fraction.num), den);

    int64_t numerator = (num / gcd1) * (fraction.num / gcd2);
    int64_t denominator = (den / gcd2) * (fraction.den / gcd1);
    return Fraction(numerator, denominator);
}

Fraction Fraction::operator/(const Fraction &fraction) const {
    if (fraction.num == 0)
        throw std::runtime_error("Fraction: division by zero");

    int64_t gcd1 = std::gcd(abs(num), fraction.num);
    int64_t gcd2 = std::gcd(abs(fraction.den), den);

    int64_t numerator = (num / gcd1) * (fraction.den / gcd2);
    int64_t denominator = (den / gcd2) * (fraction.num / gcd1);
    return Fraction(numerator, denominator);
}

Fraction& Fraction::operator+=(const Fraction &fraction) {
    *this = *this + fraction;
    return *this;
}

Fraction& Fraction::operator-=(const Fraction &fraction) {
    *this = *this - fraction;
    return *this;
}

Fraction& Fraction::operator*=(const Fraction &fraction) {
    *this = *this * fraction;
    return *this;
}

Fraction& Fraction::operator/=(const Fraction &fraction) {
    *this = *this / fraction;
    return *this;
}

bool Fraction::operator==(const Fraction &fraction) const {
    return num == fraction.num && den == fraction.den;
}

bool Fraction::operator==(int value) const {
    return num == value && den == 1;
}

bool Fraction::operator!=(const Fraction &fraction) const {
    return num != fraction.num || den != fraction.den;
}

bool Fraction::operator!=(int value) const {
    return num != value || den != 1;
}

void Fraction::normalize() {
    if (num == 0) {
        den = 1;
        return;
    }

    int64_t gcd = std::gcd(abs(num), den);

    if (gcd > 1) {
        num /= gcd;
        den /= gcd;
    }
}

std::ostream& operator<<(std::ostream &os, const Fraction &fraction) {
    return os << fraction.num << "/" << fraction.den;
}
