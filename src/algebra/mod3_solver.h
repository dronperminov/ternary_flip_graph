#pragma once

#include <iostream>
#include <vector>
#include <cstdint>

#include "../entities/mod3_vector.hpp"

class Mod3Solver {
    uint64_t rows;
    uint64_t columns;
    std::vector<uint8_t> values;
public:
    Mod3Solver(uint64_t rows, uint64_t columns);

    void set(int row, int column, uint8_t value);

    bool solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x);
};
