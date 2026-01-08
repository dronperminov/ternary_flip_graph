#pragma once

#include <iostream>
#include <vector>
#include <random>

class BinaryMatrix {
    int rows;
    int columns;
    std::vector<uint8_t> values;
public:
    BinaryMatrix(int rows, int columns);

    uint8_t operator()(int i, int j) const;
    uint8_t& operator()(int i, int j);
    uint8_t operator[](int index) const;
    uint8_t& operator[](int index);

    bool invertible(BinaryMatrix &inverse) const;

    void swapRows(int row1, int row2);
    void sandwich(const BinaryMatrix &left, const BinaryMatrix &right);
    void random(std::mt19937 &generator);
    void randomInvertible(BinaryMatrix &inverse, std::mt19937 &generator);
};
