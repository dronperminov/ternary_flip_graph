#pragma once

#include <iostream>
#include <vector>
#include <cassert>

#include "../entities/mod3_vector.hpp"

class Mod3Solver {
    int rows;
    int columns;
    std::vector<uint8_t> values;
public:
    Mod3Solver(int rows, int columns);

    void set(int row, int column, uint8_t value);

    bool solve(const std::vector<uint8_t> &b, std::vector<uint8_t> &x);
};
