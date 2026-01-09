#include "binary_lifter.h"

BinaryLifter::BinaryLifter(int n1, int n2, int n3, int rank, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, const BinarySolver &J) : jakobian(J) {
    this->dimension[0] = n1;
    this->dimension[1] = n2;
    this->dimension[2] = n3;

    this->elements[0] = n1 * n2;
    this->elements[1] = n2 * n3;
    this->elements[2] = n3 * n1;
    this->rank = rank;

    this->u = u;
    this->v = v;
    this->w = w;

    mod = 2;
    exponent = 1;
    bound = 1;

    initTensors();
}

bool BinaryLifter::lift() {
    if (exponent > 1)
        evaluateTensor();

    for (int i = 0; i < tensorSize; i++)
        b[i] = ((T0[i] - E[i]) >> exponent) & 1;

    if (!jakobian.solve(b, x))
        return false;

    updateFactor(u, elements[0], x, 0);
    updateFactor(v, elements[1], x, elements[0] * rank);
    updateFactor(w, elements[2], x, (elements[0] + elements[1]) * rank);

    exponent++;
    mod *= 2;
    bound = std::sqrt(mod / 2.0);
    return true;
}

bool BinaryLifter::reconstruct(FractionalScheme &lifted) {
    return lifted.reconstruct(dimension[0], dimension[1], dimension[2], rank, u, v, w, mod, bound);
}

int64_t BinaryLifter::getMod() const {
    return mod;
}

int64_t BinaryLifter::getBound() const {
    return bound;
}

int BinaryLifter::getExponent() const {
    return exponent;
}

void BinaryLifter::initTensors() {
    tensorSize = elements[0] * elements[1] * elements[2];
    variables = rank * (elements[0] + elements[1] + elements[2]);

    T0.resize(tensorSize);
    E.resize(tensorSize);
    evaluateTensor();

    for (int i = 0; i < tensorSize; i++)
        T0[i] = E[i] & 1;

    b.resize(tensorSize);
    x.resize(variables);
}

void BinaryLifter::evaluateTensor() {
    E.assign(tensorSize, 0);

    for (int index = 0; index < rank; index++)
        for (int i = 0; i < elements[0]; i++)
            for (int j = 0; j < elements[1]; j++)
                for (int k = 0; k < elements[2]; k++)
                    E[(i * elements[1] + j) * elements[2] + k] += u[index * elements[0] + i] * v[index * elements[1] + j] * w[index * elements[2] + k];
}

void BinaryLifter::updateFactor(std::vector<uint64_t> &f, int size, const std::vector<uint8_t> &x, int offset) {
    uint64_t mask = (uint64_t(1) << (exponent + 1)) - 1;

    for (int i = 0; i < size; i++) {
        for (int index = 0; index < rank; index++) {
            f[index * size + i] += uint64_t(x[offset + i * rank + index]) << exponent;
            f[index * size + i] &= mask;
        }
    }
}
