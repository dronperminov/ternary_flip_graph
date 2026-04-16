#include "binary_scheme.hpp"
#include "../entities/uint256_t.h"

template <typename T>
int BinaryScheme<T>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += __builtin_popcountll(uvw[i][index]);

    return count - 2 * rank - elements[2];
}

template <>
int BinaryScheme<__uint128_t>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += __builtin_popcountll(uvw[i][index]) + __builtin_popcountll(uvw[i][index] >> 64);

    return count - 2 * rank - elements[2];
}

template <>
int BinaryScheme<uint256_t>::getComplexity() const {
    int count = 0;

    for (int i = 0; i < 3; i++)
        for (int index = 0; index < rank; index++)
            count += uvw[i][index].popcount();

    return count - 2 * rank - elements[2];
}

template class BinaryScheme<uint16_t>;
template class BinaryScheme<uint32_t>;
template class BinaryScheme<uint64_t>;
template class BinaryScheme<__uint128_t>;
template class BinaryScheme<uint256_t>;
