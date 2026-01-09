#include "fraction.h"

Fraction abs(const Fraction &fraction) {
    return Fraction(abs(fraction.numerator()), fraction.denominator());
}

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

std::string Fraction::pretty() const {
    std::stringstream ss;
    ss << num;

    if (den > 1)
        ss << "/" << den;

    return ss.str();
}

bool Fraction::isInteger() const {
    return den == 1;
}

bool Fraction::isTernaryInteger() const {
    return den == 1 && -1 <= num && num <= 1;
}

bool Fraction::reconstruct(int64_t a, int64_t mod, int64_t bound) {
    a = ((a % mod) + mod) % mod;

    int64_t r0 = mod;
    int64_t r1 = a;
    int64_t t0 = 0;
    int64_t t1 = 1;

    while (r1 != 0 && r1 > bound) {
        int64_t q = r0 / r1;
        int64_t r2 = r0 - q * r1;
        int64_t t2 = t0 - q * t1;

        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }

    if (abs(r1) > bound || abs(t1) > bound || t1 == 0)
        return false;

    if (t1 < 0) {
        r1 = -r1;
        t1 = -t1;
    }

    if (std::gcd(abs(r1), abs(t1)) != 1)
        return false;

    num = r1;
    den = t1;
    normalize();
    return true;
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

bool Fraction::operator>(const Fraction& fraction) const {
    return num * fraction.den > fraction.num * den;
}

bool Fraction::operator<(const Fraction& fraction) const {
    return num * fraction.den < fraction.num * den;
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

std::istream& operator>>(std::istream &is, Fraction &fraction) {
    is >> fraction.num >> fraction.den;
    fraction.normalize();
    return is;
}
