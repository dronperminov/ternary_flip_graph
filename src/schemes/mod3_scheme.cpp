#include "mod3_scheme.hpp"
#include "../entities/uint256_t.h"

template class Mod3Scheme<uint16_t>;
template class Mod3Scheme<uint32_t>;
template class Mod3Scheme<uint64_t>;
template class Mod3Scheme<__uint128_t>;
template class Mod3Scheme<uint256_t>;
