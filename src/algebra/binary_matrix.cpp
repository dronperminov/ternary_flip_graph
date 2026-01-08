#include "binary_matrix.h"

BinaryMatrix::BinaryMatrix(int rows, int columns) {
    this->rows = rows;
    this->columns = columns;
    this->values.assign(rows * columns, 0);
}

uint8_t BinaryMatrix::operator()(int i, int j) const {
    return values[i * columns + j];
}

uint8_t& BinaryMatrix::operator()(int i, int j) {
    return values[i * columns + j];
}

uint8_t BinaryMatrix::operator[](int index) const {
    return values[index];
}

uint8_t& BinaryMatrix::operator[](int index) {
    return values[index];
}

bool BinaryMatrix::invertible(BinaryMatrix &inverse) const {
    if (rows != columns)
        return false;

    int size = rows;
    int size2 = size * 2;
    BinaryMatrix augmented(size, size2);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            augmented(i, j) = values[i * size + j];
            augmented(i, j + size) = i == j;
        }
    }

    for (int col = 0; col < size; col++) {
        int pivotRow = col;

        while (pivotRow < size && !augmented(pivotRow, col))
            pivotRow++;

        if (pivotRow == size)
            return false;

        if (pivotRow != col)
            augmented.swapRows(col, pivotRow);

        for (int i = 0; i < size; i++)
            if (i != col && augmented(i, col))
                for (int j = 0; j < size2; j++)
                    augmented(i, j) ^= augmented(col, j);
    }

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            inverse(i, j) = augmented(i, j + size);

    return true;
}

void BinaryMatrix::swapRows(int row1, int row2) {
    for (int i = 0; i < columns; i++)
        std::swap(values[row1 * columns + i], values[row2 * columns + i]);
}

void BinaryMatrix::sandwich(const BinaryMatrix &left, const BinaryMatrix &right) {
    BinaryMatrix tmp(rows, columns);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            uint8_t value = 0;

            for (int k = 0; k < rows; k++)
                value ^= left.values[i * rows + k] & values[k * columns + j];

            tmp.values[i * columns + j] = value;
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            uint8_t value = 0;

            for (int k = 0; k < columns; k++)
                value ^= tmp.values[i * columns + k] & right.values[k * columns + j];

            values[i * columns + j] = value;
        }
    }
}

void BinaryMatrix::random(std::mt19937 &generator) {
    for (int i = 0; i < rows * columns; i++)
        values[i] = generator() & 1;
}

void BinaryMatrix::randomInvertible(BinaryMatrix &inverse, std::mt19937 &generator) {
    do {
        random(generator);
    } while (!invertible(inverse));
}
