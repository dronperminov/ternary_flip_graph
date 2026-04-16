#include "mod3_vector.hpp"
#include "uint256_t.h"

template <typename T>
int Mod3Vector<T>::nonZeroCount() const {
    return __builtin_popcountll(low | high);
}

template <>
int Mod3Vector<__uint128_t>::nonZeroCount() const {
    __uint128_t values = low | high;
    return __builtin_popcountll(values) + __builtin_popcountll(values >> 64);
}

template <>
int Mod3Vector<uint256_t>::nonZeroCount() const {
    uint256_t values = low | high;
    return values.popcount();
}

template class Mod3Vector<uint16_t>;
template class Mod3Vector<uint32_t>;
template class Mod3Vector<uint64_t>;
template class Mod3Vector<__uint128_t>;
template class Mod3Vector<uint256_t>;
