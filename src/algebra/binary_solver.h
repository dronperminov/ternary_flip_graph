#pragma once

#include <iostream>
#include <vector>

class BinarySolver {
    int rows;
    int columns;
    std::vector<uint8_t> values;
public:
    BinarySolver(int rows, int columns);

    void inverse(int row, int column);

    bool solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x);
};
