#include "mod3_solver.h"

Mod3Solver::Mod3Solver(uint64_t rows, uint64_t columns) : wordsPerRow(columns / 64 + 1), values(rows * wordsPerRow, 0) {
    this->rows = rows;
    this->columns = columns;
}

void Mod3Solver::set(int row, int column, uint8_t value) {
    values[row * wordsPerRow + column / 64].set(column % 64, value % 3);
}

bool Mod3Solver::solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x) {
    std::vector<Mod3Vector<uint64_t>> augmented(values);
    std::vector<int> pivotCol(rows, -1);

    uint64_t lastWord = columns / 64;
    uint64_t lastBit = columns % 64;

    for (uint64_t i = 0; i < rows; i++)
        augmented[i * wordsPerRow + lastWord].set(lastBit, b[i]);

    uint64_t rank = 0;

    for (uint64_t column = 0; column < columns && rank < rows; column++) {
        uint64_t word = column / 64;
        uint64_t bit = column % 64;

        uint64_t pivotRow = rank;
        while (pivotRow < rows && !augmented[pivotRow * wordsPerRow + word][bit])
            pivotRow++;

        if (pivotRow == rows)
            continue;

        if (pivotRow != rank) {
            uint64_t offset1 = rank * wordsPerRow + word;
            uint64_t offset2 = pivotRow * wordsPerRow + word;

            for (uint64_t j = word; j < wordsPerRow; j++)
                std::swap(augmented[offset1++], augmented[offset2++]);
        }

        int pivotValue = augmented[rank * wordsPerRow + word][bit];
        if (pivotValue != 1) {
            uint64_t offset = rank * wordsPerRow + word;

            for (uint64_t j = word; j < wordsPerRow; j++)
                augmented[offset++] *= pivotValue;
        }

        for (uint64_t row = 0; row < rows; row++) {
            int pivot = augmented[row * wordsPerRow + word][bit];

            if (row != rank && pivot) {
                uint64_t offset1 = row * wordsPerRow + word;
                uint64_t offset2 = rank * wordsPerRow + word;

                for (uint64_t j = word; j < wordsPerRow; j++)
                    augmented[offset1++] -= augmented[offset2++] * pivot;
            }
        }

        pivotCol[rank++] = column;
    }

    for (uint64_t i = rank; i < rows; i++)
        if (augmented[i * wordsPerRow + lastWord][lastBit])
            return false;

    x.assign(columns, 0);
    for (uint64_t i = 0; i < rank; i++)
        x[pivotCol[i]] = augmented[i * wordsPerRow + lastWord][lastBit];

    return true;
}
