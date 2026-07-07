#pragma once

#include <iostream>
#include <cstdint>
#include <vector>
#include "utils.h"

class ModMatrix {
    int rows;
    int columns;
    int64_t mod;
    std::vector<int64_t> values;
public:
    ModMatrix(int rows, int columns, int64_t mod);

    int64_t& operator()(int i, int j);
    int64_t operator()(int i, int j) const;

    int64_t& operator[](int index);
    int64_t operator[](int index) const;

    int rank() const;

    void swapRows(int row1, int row2, int column = 0);
    void multiplyRow(int row, int64_t multiplier, int column = 0);
    void subtractRow(int row1, int row2, int64_t value, int column = 0);
};
