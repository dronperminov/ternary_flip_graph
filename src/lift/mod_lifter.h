#pragma once

#include <iostream>
#include <vector>

#include "../algebra/mod_matrix.h"
#include "../schemes/fractional_scheme.h"

class ModLifter {
    int dimension[3];
    int elements[3];
    int rank;
    std::vector<uint64_t> u;
    std::vector<uint64_t> v;
    std::vector<uint64_t> w;

    int tensorSize;
    int variables;
    int64_t p;
    int64_t mod;
    int64_t bound;
    int exponent;

    std::vector<int64_t> T0;
    std::vector<int64_t> E;

    ModMatrix jakobian;
    std::vector<int64_t> b;
    std::vector<int64_t> x;
public:
    ModLifter(int n1, int n2, int n3, int rank, int64_t p, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, const ModMatrix &jakobian);

    bool lift();
    bool canLift();
    bool reconstruct(FractionalScheme &lifted);

    int64_t getMod() const;
    int64_t getBound() const;
    int getExponent() const;

    void show() const;
private:
    void initTensors();
    void evaluateTensor();
    void updateFactor(std::vector<uint64_t> &f, int size, const std::vector<int64_t> &x, int offset);
};
