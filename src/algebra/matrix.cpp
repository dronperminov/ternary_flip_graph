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

    for (int column = 0; column < size; column++) {
        int pivotRow = column;

        for (int row = column + 1; row < size; row++)
            if (abs(augmented(row, column)) > abs(augmented(pivotRow, column)))
                pivotRow = row;

        if (pivotRow != column)
            augmented.swapRows(column, pivotRow);

        if (augmented(column, column) == 0)
            return false;

        Fraction pivot = augmented(column, column);

        for (int j = 0; j < size * 2; j++)
            augmented(column, j) /= pivot;

        for (int i = 0; i < size; i++) {
            if (i == column)
                continue;

            Fraction factor = augmented(i, column);

            for (int j = 0; j < size * 2; j++)
                augmented(i, j) -= factor * augmented(column, j);
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

int Matrix::rank() const {
    Matrix tmp(rows, columns);
    tmp.values = values;

    int rank = 0;

    for (int column = 0; column < columns && rank < rows; column++) {
        int pivotRow = rank;
        while (pivotRow < rows && tmp(pivotRow, column) == 0)
            pivotRow++;

        if (pivotRow == rows)
            continue;

        if (pivotRow != rank)
            tmp.swapRows(pivotRow, rank, column);

        Fraction pivot = tmp(rank, column);
        tmp.divideRow(rank, pivot, column);

        for (int row = rank + 1; row < rows; row++) {
            Fraction value = tmp(row, column);
            tmp.subtractRow(row, rank, value, column);
        }

        rank++;
    }

    return rank;
}

void Matrix::swapRows(int row1, int row2, int column) {
    for (int i = column; i < columns; i++)
        std::swap(values[row1 * columns + i], values[row2 * columns + i]);
}

void Matrix::divideRow(int row, const Fraction &divider, int column) {
    for (int i = column; i < columns; i++)
        values[row * columns + i] /= divider;
}

void Matrix::subtractRow(int row1, int row2, const Fraction &value, int column) {
    for (int i = column; i < columns; i++)
        values[row1 * columns + i] -= values[row2 * columns + i] * value;
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
            os << matrix(i, j).pretty() << " ";

        os << std::endl;
    }

    return os;
}
