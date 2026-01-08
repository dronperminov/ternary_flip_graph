#pragma once

#include <iostream>

class BinarySolver {
    int rows;
    int columns;
    std::vector<uint8_t> values;
public:
    BinarySolver(int rows, int columns);

    void inverse(int row, int column);

    bool solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x);
};

BinarySolver::BinarySolver(int rows, int columns) {
    this->rows = rows;
    this->columns = columns;

    values.assign(rows * columns, 0);
}

void BinarySolver::inverse(int row, int column) {
    values[row * columns + column] ^= 1;
}

bool BinarySolver::solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x) {
    int augmentedColumns = columns + 1;
    std::vector<uint8_t> augmented(rows * augmentedColumns, 0);
    std::vector<int> pivotCol(rows, -1);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++)
            augmented[i * augmentedColumns + j] = values[i * columns + j];

        augmented[i * augmentedColumns + columns] = b[i] ? 1 : 0;
    }

    int rank = 0;

    for (int column = 0; column < columns && rank < rows; column++) {
        int pivotRow = -1;
        for (int row = rank; row < rows; row++) {
            if (augmented[row * augmentedColumns + column]) {
                pivotRow = row;
                break;
            }
        }

        if (pivotRow == -1)
            continue;

        if (pivotRow != rank) {
            int offset1 = rank * augmentedColumns + column;
            int offset2 = pivotRow * augmentedColumns + column;

            for (int j = column; j < augmentedColumns; j++)
                std::swap(augmented[offset1++], augmented[offset2++]);
        }

        for (int row = 0; row < rows; row++) {
            if (row != rank && augmented[row * augmentedColumns + column]) {
                int offset1 = row * augmentedColumns + column;
                int offset2 = rank * augmentedColumns + column;

                for (int j = column; j < augmentedColumns; j++)
                    augmented[offset1++] ^= augmented[offset2++];
            }
        }

        pivotCol[rank++] = column;
    }

    for (int i = rank; i < rows; i++)
        if (augmented[i * augmentedColumns + columns] != 0)
            return false;

    x.assign(columns, 0);
    for (int i = 0; i < rank; i++)
        x[pivotCol[i]] = augmented[i * augmentedColumns + columns];

    return true;
}
