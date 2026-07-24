#include "mod_lifter.h"

ModLifter::ModLifter(int n1, int n2, int n3, int rank, int64_t p, const std::vector<uint64_t> &u, const std::vector<uint64_t> &v, const std::vector<uint64_t> &w, const ModMatrix &J) : jakobian(J) {
    this->dimension[0] = n1;
    this->dimension[1] = n2;
    this->dimension[2] = n3;

    this->elements[0] = n1 * n2;
    this->elements[1] = n2 * n3;
    this->elements[2] = n3 * n1;
    this->rank = rank;
    this->p = p;

    this->u = u;
    this->v = v;
    this->w = w;

    mod = p;
    exponent = 1;
    bound = 1;

    initTensors();
}

bool ModLifter::lift() {
    if (exponent > 1)
        evaluateTensor();

    for (int i = 0; i < tensorSize; i++)
        b[i] = (((T0[i] - E[i]) / mod) % p + p) % p;

    if (!jakobian.solve(b, x))
        return false;

    updateFactor(u, elements[0], x, 0);
    updateFactor(v, elements[1], x, elements[0] * rank);
    updateFactor(w, elements[2], x, (elements[0] + elements[1]) * rank);

    exponent++;
    mod *= p;
    bound = std::sqrt(mod / 2.0);
    return true;
}

bool ModLifter::canLift() {
    for (int i = 0; i < tensorSize; i++)
        b[i] = (((T0[i] - E[i]) / mod) % p + p) % p;

    return jakobian.solve(b, x);
}

bool ModLifter::reconstruct(FractionalScheme &lifted) {
    return lifted.reconstruct(dimension[0], dimension[1], dimension[2], rank, u, v, w, mod, bound);
}

int64_t ModLifter::getMod() const {
    return mod;
}

int64_t ModLifter::getBound() const {
    return bound;
}

int ModLifter::getExponent() const {
    return exponent;
}

void ModLifter::show() const {
    std::cout << "mod: " << mod;
    std::cout << std::endl << "U:";
    for (int i = 0; i < rank * elements[0]; i++)
        std::cout << " " << u[i];

    std::cout << std::endl << "V:";
    for (int i = 0; i < rank * elements[1]; i++)
        std::cout << " " << v[i];

    std::cout << std::endl << "W:";
    for (int i = 0; i < rank * elements[2]; i++)
        std::cout << " " << w[i];

    std::cout << std::endl;
}

void ModLifter::initTensors() {
    tensorSize = elements[0] * elements[1] * elements[2];
    variables = rank * (elements[0] + elements[1] + elements[2]);

    T0.resize(tensorSize);
    E.resize(tensorSize);
    evaluateTensor();

    for (int i = 0; i < tensorSize; i++)
        T0[i] = E[i] % p;

    b.resize(tensorSize);
    x.resize(variables);
}

void ModLifter::evaluateTensor() {
    E.assign(tensorSize, 0);

    for (int index = 0; index < rank; index++)
        for (int i = 0; i < elements[0]; i++)
            for (int j = 0; j < elements[1]; j++)
                for (int k = 0; k < elements[2]; k++)
                    E[(i * elements[1] + j) * elements[2] + k] += u[index * elements[0] + i] * v[index * elements[1] + j] * w[index * elements[2] + k];
}

void ModLifter::updateFactor(std::vector<uint64_t> &f, int size, const std::vector<int64_t> &x, int offset) {
    uint64_t modNext = mod * p;

    for (int i = 0; i < size; i++)
        for (int index = 0; index < rank; index++)
            f[index * size + i] = (f[index * size + i] + x[offset + i * rank + index] * mod) % modNext;
}
