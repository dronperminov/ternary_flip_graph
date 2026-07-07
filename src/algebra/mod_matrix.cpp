#include "mod_matrix.h"

ModMatrix::ModMatrix(int rows, int columns, int64_t mod) {
    this->rows = rows;
    this->columns = columns;
    this->mod = mod;
    this->values.resize(rows * columns);
}

int64_t& ModMatrix::operator()(int i, int j) {
    return values[i * columns + j];
}

int64_t ModMatrix::operator()(int i, int j) const {
    return values[i * columns + j];
}

int64_t& ModMatrix::operator[](int index) {
    return values[index];
}

int64_t ModMatrix::operator[](int index) const {
    return values[index];
}

int ModMatrix::rank() const {
    ModMatrix tmp(rows, columns, mod);
    tmp.values = values;

    int rank = 0;
    std::vector<int> rowIndices(rows);
    for (int i = 0; i < rows; i++)
        rowIndices[i] = i;

    for (int column = 0; column < columns && rank < rows; column++) {
        int pivotRow = rank;
        while (pivotRow < rows && tmp(rowIndices[pivotRow], column) == 0)
            pivotRow++;

        if (pivotRow == rows)
            continue;

        if (pivotRow != rank)
            std::swap(rowIndices[pivotRow], rowIndices[rank]);

        int64_t pivot = modInverse(tmp(rowIndices[rank], column), mod);
        tmp.multiplyRow(rowIndices[rank], pivot, column);

        for (int row = rank + 1; row < rows; row++) {
            int64_t value = tmp(rowIndices[row], column);
            tmp.subtractRow(rowIndices[row], rowIndices[rank], value, column);
        }

        rank++;
    }

    return rank;
}

void ModMatrix::swapRows(int row1, int row2, int column) {
    for (int i = column; i < columns; i++)
        std::swap(values[row1 * columns + i], values[row2 * columns + i]);
}

void ModMatrix::multiplyRow(int row, int64_t multiplier, int column) {
    for (int i = column; i < columns; i++)
        values[row * columns + i] = (values[row * columns + i] * multiplier) % mod;
}

void ModMatrix::subtractRow(int row1, int row2, int64_t value, int column) {
    for (int i = column; i < columns; i++)
        values[row1 * columns + i] = (values[row1 * columns + i] + mod - (values[row2 * columns + i] * value) % mod) % mod;
}
