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
    int wordsPerRow = (columns + 1 + 63) / 64;
    std::vector<uint64_t> augmented(rows * wordsPerRow, 0);
    std::vector<int> pivotCol(rows, -1);

    int lastWord = columns / 64;
    uint64_t lastMask = uint64_t(1) << (columns % 64);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++)
            augmented[i * wordsPerRow + j / 64] |= uint64_t(values[i * columns + j]) << (j % 64);

        if (b[i])
            augmented[i * wordsPerRow + lastWord] |= lastMask;
    }

    int rank = 0;

    for (int column = 0; column < columns && rank < rows; column++) {
        int word = column / 64;
        uint64_t mask = uint64_t(1) << (column % 64);

        int pivotRow = -1;
        for (int row = rank; row < rows; row++) {
            if (augmented[row * wordsPerRow + word] & mask) {
                pivotRow = row;
                break;
            }
        }

        if (pivotRow == -1)
            continue;

        if (pivotRow != rank) {
            int offset1 = rank * wordsPerRow + word;
            int offset2 = pivotRow * wordsPerRow + word;

            for (int j = word; j < wordsPerRow; j++)
                std::swap(augmented[offset1++], augmented[offset2++]);
        }

        for (int row = 0; row < rows; row++) {
            if (row != rank && (augmented[row * wordsPerRow + word] & mask)) {
                int offset1 = row * wordsPerRow + word;
                int offset2 = rank * wordsPerRow + word;

                for (int j = word; j < wordsPerRow; j++)
                    augmented[offset1++] ^= augmented[offset2++];
            }
        }

        pivotCol[rank++] = column;
    }

    for (int i = rank; i < rows; i++)
        if (augmented[i * wordsPerRow + lastWord] & lastMask)
            return false;

    x.assign(columns, 0);
    for (int i = 0; i < rank; i++)
        x[pivotCol[i]] = augmented[i * wordsPerRow + lastWord] & lastMask ? 1 : 0;

    return true;
}
