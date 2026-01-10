#pragma once

#include <iostream>
#include <vector>

class BinarySolver {
    int rows;
    int columns;
    std::vector<uint8_t> values;
    std::vector<int8_t> xs;
public:
    BinarySolver(int rows, int columns);

    void set(int row, int column, uint8_t value);
    void setVariable(int variable, uint8_t value);
    void reset();

    bool solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x);
};
