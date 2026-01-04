#include "matrix.h"

Matrix::Matrix(int rows, int columns) {
    this->rows = rows;
    this->columns = columns;
    this->values.resize(rows * columns);
}

Fraction& Matrix::operator()(int i, int j) {
    return values[i * columns + j];
}

Fraction Matrix::operator()(int i, int j) const {
    return values[i * columns + j];
}

Fraction& Matrix::operator[](int index) {
    return values[index];
}

Fraction Matrix::operator[](int index) const {
    return values[index];
}

Matrix Matrix::operator*(const Matrix &matrix) const {
    if (columns != matrix.rows)
        throw std::runtime_error("matrix sizes mismatch");

    Matrix result(rows, matrix.columns);

    for (int i = 0; i < result.rows; i++) {
        for (int j = 0; j < result.columns; j++) {
            Fraction sum = 0;

            for (int k = 0; k < columns; k++)
                sum += values[i * columns + k] * matrix.values[k * matrix.columns + j];

            result.values[i * result.columns + j] = sum;
        }
    }

    return result;
}

bool Matrix::invertible(Matrix &inverse) const {
    if (rows != columns)
        return false;

    int size = rows;
    Matrix augmented(size, size * 2);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            augmented(i, j) = values[i * size + j];
            augmented(i, j + size) = i == j ? 1 : 0;
        }
    }

    for (int col = 0; col < size; col++) {
        int pivotRow = col;

        for (int row = col + 1; row < size; row++)
            if (abs(augmented(row, col)) > abs(augmented(pivotRow, col)))
                pivotRow = row;

        if (pivotRow != col)
            augmented.swapRows(col, pivotRow);

        if (augmented(col, col) == 0)
            return false;

        Fraction pivot = augmented(col, col);

        for (int j = 0; j < size * 2; j++)
            augmented(col, j) /= pivot;

        for (int i = 0; i < size; i++) {
            if (i == col)
                continue;

            Fraction factor = augmented(i, col);

            for (int j = 0; j < size * 2; j++)
                augmented(i, j) -= factor * augmented(col, j);
        }
    }

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            inverse(i, j) = augmented(i, j + size);

    return true;
}

bool Matrix::isTernary() const {
    for (int i = 0; i < rows * columns; i++)
        if (!values[i].isTernaryInteger())
            return false;

    return true;
}

int Matrix::fractionsCount() const {
    int count = 0;

    for (int i = 0; i < rows * columns; i++)
        if (!values[i].isInteger())
            count++;

    return count;
}

void Matrix::swapRows(int row1, int row2) {
    for (int i = 0; i < columns; i++)
        std::swap(values[row1 * columns + i], values[row2 * columns + i]);
}

void Matrix::sandwich(const Matrix &left, const Matrix &right) {
    Matrix tmp(rows, columns);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            Fraction value = 0;

            for (int k = 0; k < rows; k++)
                value += left.values[i * rows + k] * values[k * columns + j];

            tmp.values[i * columns + j] = value;
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            Fraction value = 0;

            for (int k = 0; k < columns; k++)
                value += tmp.values[i * columns + k] * right.values[k * columns + j];

            values[i * columns + j] = value;
        }
    }
}

void Matrix::random(int min, int max, int denominator, std::mt19937 &generator) {
    std::uniform_int_distribution<int> distribution(min, max);

    for (int i = 0; i < rows * columns; i++)
        values[i] = Fraction(distribution(generator), denominator);
}

void Matrix::diagonal(const Fraction &value) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++)
            values[i * columns + j] = 0;

        values[i * columns + i] = value;
    }
}

bool Matrix::toRing(int ring) {
    for (int i = 0; i < rows * columns; i++) {
        int a = ((values[i].numerator() % ring) + ring) % ring;
        int b = values[i].denominator();
        int c = 0;

        while (c < ring && ((b * c) % ring != a))
            c++;

        if (c == ring)
            return false;

        values[i] = c;
    }

    return true;
}

std::istream& operator>>(std::istream& is, Matrix &matrix) {
    for (int i = 0; i < matrix.rows * matrix.columns; i++)
        is >> matrix.values[i];

    return is;
}

std::ostream& operator<<(std::ostream& os, const Matrix &matrix) {
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.columns; j++)
            os << matrix(i, j) << " ";

        os << std::endl;
    }

    return os;
}
