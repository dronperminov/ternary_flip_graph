#pragma once

#include <iostream>
#include <vector>

class BinarySolver {
    uint64_t rows;
    uint64_t columns;
    uint64_t wordsPerRow;
    std::vector<uint64_t> values;
    std::vector<int8_t> xs;
public:
    BinarySolver(uint64_t rows, uint64_t columns);

    void set(int row, int column, uint8_t value);
    void setVariable(int variable, uint8_t value);
    void reset();

    bool solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x);
};
