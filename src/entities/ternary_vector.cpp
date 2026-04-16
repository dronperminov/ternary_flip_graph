#include "ternary_vector.hpp"
#include "uint256_t.h"

template <typename T>
int TernaryVector<T>::nonZeroCount() const {
    return __builtin_popcountll(values);
}

template <>
int TernaryVector<__uint128_t>::nonZeroCount() const {
    return __builtin_popcountll(values) + __builtin_popcountll(values >> 64);
}

template <>
int TernaryVector<uint256_t>::nonZeroCount() const {
    return values.popcount();
}

template class TernaryVector<uint16_t>;
template class TernaryVector<uint32_t>;
template class TernaryVector<uint64_t>;
template class TernaryVector<__uint128_t>;
template class TernaryVector<uint256_t>;
