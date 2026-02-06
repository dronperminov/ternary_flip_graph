#include "binary_solver.h"
#include <cassert>

BinarySolver::BinarySolver(uint64_t rows, uint64_t columns) : wordsPerRow(columns / 64 + 1), values(rows * wordsPerRow, 0), xs(columns, -1) {
    this->rows = rows;
    this->columns = columns;
}

void BinarySolver::set(int row, int column, uint8_t value) {
    uint64_t mask = uint64_t(1) << (column % 64);

    if (value & 1)
        values[row * wordsPerRow + column / 64] |= mask;
    else
        values[row * wordsPerRow + column / 64] &= ~mask;
}

void BinarySolver::setVariable(int variable, uint8_t value) {
    xs[variable] = value;
}

void BinarySolver::reset() {
    xs.assign(columns, -1);
}

bool BinarySolver::solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x) {
    std::vector<uint64_t> augmented(rows * wordsPerRow, 0);
    std::vector<int> pivotCol(rows, -1);

    uint64_t lastWord = columns / 64;
    uint64_t lastMask = uint64_t(1) << (columns % 64);

    for (uint64_t i = 0; i < rows; i++) {
        uint8_t bi = b[i];

        for (uint64_t j = 0; j < columns; j++) {
            uint8_t xj = (values[i * wordsPerRow + j / 64] >> (j % 64)) & 1;
            augmented[i * wordsPerRow + j / 64] |= uint64_t(xj) << (j % 64);

            if (xj && xs[j] > -1)
                bi ^= xs[j];
        }

        if (bi)
            augmented[i * wordsPerRow + lastWord] |= lastMask;
    }

    uint64_t rank = 0;

    x.assign(columns, 0);

    for (uint64_t column = 0; column < columns && rank < rows; column++) {
        if (xs[column] > -1) {
            x[column] = xs[column];
            continue;
        }

        uint64_t word = column / 64;
        uint64_t mask = uint64_t(1) << (column % 64);

        uint64_t pivotRow = rank;
        while (pivotRow < rows && !(augmented[pivotRow * wordsPerRow + word] & mask))
            pivotRow++;

        if (pivotRow == rows)
            continue;

        if (pivotRow != rank) {
            uint64_t offset1 = rank * wordsPerRow + word;
            uint64_t offset2 = pivotRow * wordsPerRow + word;

            for (uint64_t j = word; j < wordsPerRow; j++)
                std::swap(augmented[offset1++], augmented[offset2++]);
        }

        for (uint64_t row = 0; row < rows; row++) {
            if (row != rank && (augmented[row * wordsPerRow + word] & mask)) {
                uint64_t offset1 = row * wordsPerRow + word;
                uint64_t offset2 = rank * wordsPerRow + word;

                for (uint64_t j = word; j < wordsPerRow; j++)
                    augmented[offset1++] ^= augmented[offset2++];
            }
        }

        pivotCol[rank++] = column;
    }

    for (uint64_t i = rank; i < rows; i++)
        if (augmented[i * wordsPerRow + lastWord] & lastMask)
            return false;

    for (uint64_t i = 0; i < rank; i++)
        x[pivotCol[i]] = augmented[i * wordsPerRow + lastWord] & lastMask ? 1 : 0;

    return true;
}
