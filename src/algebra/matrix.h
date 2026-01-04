#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <cassert>

#include "fraction.h"

class Matrix {
    int rows;
    int columns;
    std::vector<Fraction> values;
public:
    Matrix(int rows, int columns);

    Fraction& operator()(int i, int j);
    Fraction operator()(int i, int j) const;

    Fraction& operator[](int index);
    Fraction operator[](int index) const;

    Matrix operator*(const Matrix &matrix) const;

    bool invertible(Matrix &inverse) const;
    bool isTernary() const;

    int fractionsCount() const;

    void swapRows(int row1, int row2);
    void sandwich(const Matrix &left, const Matrix &right);
    void random(int min, int max, int denominator, std::mt19937 &generator);
    void diagonal(const Fraction &value);

    bool toRing(int ring);

    friend std::istream& operator>>(std::istream& is, Matrix &matrix);
    friend std::ostream& operator<<(std::ostream& os, const Matrix &matrix);
};
