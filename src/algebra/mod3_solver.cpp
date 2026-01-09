#include "mod3_solver.h"

Mod3Solver::Mod3Solver(int rows, int columns) {
    this->rows = rows;
    this->columns = columns;

    values.assign(rows * columns, 0);
}

void Mod3Solver::set(int row, int column, uint8_t value) {
    values[row * columns + column] = value % 3;
}

bool Mod3Solver::solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x) {
    int wordsPerRow = (columns + 1 + 63) / 64;
    std::vector<Mod3Vector<uint64_t>> augmented(rows * wordsPerRow, 64);
    std::vector<int> pivotCol(rows, -1);

    int lastWord = columns / 64;
    int lastBit = columns % 64;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++)
            augmented[i * wordsPerRow + j / 64].set(j % 64, values[i * columns + j]);

        augmented[i * wordsPerRow + lastWord].set(lastBit, b[i] % 3);
    }

    int rank = 0;

    for (int column = 0; column < columns && rank < rows; column++) {
        int word = column / 64;
        int bit = column % 64;

        int pivotRow = -1;
        for (int row = rank; row < rows; row++) {
            if (augmented[row * wordsPerRow + word][bit]) {
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

        int pivotValue = augmented[rank * wordsPerRow + word][bit];
        if (pivotValue != 1) {
            int offset = rank * wordsPerRow + word;

            for (int j = word; j < wordsPerRow; j++)
                augmented[offset++] *= pivotValue;
        }

        for (int row = 0; row < rows; row++) {
            int pivot = augmented[row * wordsPerRow + word][bit];

            if (row != rank && pivot) {
                int offset1 = row * wordsPerRow + word;
                int offset2 = rank * wordsPerRow + word;

                for (int j = word; j < wordsPerRow; j++)
                    augmented[offset1++] -= augmented[offset2++] * pivot;
            }
        }

        pivotCol[rank++] = column;
    }

    for (int i = rank; i < rows; i++)
        if (augmented[i * wordsPerRow + lastWord][lastBit])
            return false;

    x.assign(columns, 0);
    for (int i = 0; i < rank; i++)
        x[pivotCol[i]] = augmented[i * wordsPerRow + lastWord][lastBit];

    return true;
}
