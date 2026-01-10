#include "binary_solver.h"
#include <cassert>

BinarySolver::BinarySolver(int rows, int columns) : values(rows * columns, 0), xs(columns, -1) {
    this->rows = rows;
    this->columns = columns;
}

void BinarySolver::set(int row, int column, uint8_t value) {
    values[row * columns + column] = value & 1;
}

void BinarySolver::setVariable(int variable, uint8_t value) {
    xs[variable] = value;
}

void BinarySolver::reset() {
    xs.assign(columns, -1);
}

bool BinarySolver::solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x) {
    int wordsPerRow = (columns + 1 + 63) / 64;
    std::vector<uint64_t> augmented(rows * wordsPerRow, 0);
    std::vector<int> pivotCol(rows, -1);

    int lastWord = columns / 64;
    uint64_t lastMask = uint64_t(1) << (columns % 64);

    for (int i = 0; i < rows; i++) {
        uint8_t bi = b[i];

        for (int j = 0; j < columns; j++) {
            uint8_t xj = values[i * columns + j];
            augmented[i * wordsPerRow + j / 64] |= uint64_t(xj) << (j % 64);

            if (xj && xs[j] > -1)
                bi ^= xs[j];
        }

        if (bi)
            augmented[i * wordsPerRow + lastWord] |= lastMask;
    }

    int rank = 0;

    x.assign(columns, 0);

    for (int column = 0; column < columns && rank < rows; column++) {
        if (xs[column] > -1) {
            x[column] = xs[column];
            continue;
        }

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

    for (int i = 0; i < rank; i++)
        x[pivotCol[i]] = augmented[i * wordsPerRow + lastWord] & lastMask ? 1 : 0;

    return true;
}
