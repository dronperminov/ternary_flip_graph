#pragma once

#include <iostream>
#include <vector>

#include "../algebra/binary_solver.h"
#include "../schemes/fractional_scheme.h"

class BinaryLifter {
    int dimension[3];
    int elements[3];
    int rank;
    std::vector<uint64_t> u;
    std::vector<uint64_t> v;
    std::vector<uint64_t> w;

    int tensorSize;
    int variables;
    int64_t mod;
    int64_t bound;
    int exponent;

    std::vector<int64_t> T0;
    std::vector<int64_t> E;

    BinarySolver jakobian;
    std::vector<uint8_t> b;
    std::vector<uint8_t> x;
public:
    BinaryLifter(int n1, int n2, int n3, int rank, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, const BinarySolver &jakobian);

    bool lift();
    bool reconstruct(FractionalScheme &lifted);

    int64_t getMod() const;
    int64_t getBound() const;
    int getExponent() const;
private:
    void initTensors();
    void evaluateTensor();
    void updateFactor(std::vector<uint64_t> &f, int size, const std::vector<uint8_t> &x, int offset);
};
