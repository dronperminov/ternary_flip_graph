#include "uint256_t.h"

uint256_t::uint256_t() {
    high = 0;
    low = 0;
}

uint256_t::uint256_t(int value) {
    high = 0;
    low = value;
}

uint256_t::uint256_t(unsigned int value) {
    high = 0;
    low = value;
}

uint256_t::uint256_t(__uint128_t high, __uint128_t low) {
    this->high = high;
    this->low = low;
}

uint256_t uint256_t::operator<<(int n) const {
    if (n == 0)
        return *this;

    if (n >= 256)
        return uint256_t(0);

    uint256_t result;

    if (n < 128) {
        result.high = (high << n) | (low >> (128 - n));
        result.low = low << n;
    } else {
        result.high = low << (n - 128);
        result.low = 0;
    }

    return result;
}

uint256_t uint256_t::operator>>(int n) const {
    if (n == 0)
        return *this;

    if (n >= 256)
        return uint256_t(0);

    uint256_t result;

    if (n < 128) {
        result.low = (low >> n) | (high << (128 - n));
        result.high = high >> n;
    } else {
        result.low = high >> (n - 128);
        result.high = 0;
    }

    return result;
}

uint256_t uint256_t::operator&(const uint256_t &other) const {
    return uint256_t(high & other.high, low & other.low);
}

uint256_t& uint256_t::operator&=(const uint256_t &other) {
    high &= other.high;
    low &= other.low;
    return *this;
}

uint256_t uint256_t::operator&(int value) const {
    return uint256_t(0, low & value);
}

uint256_t uint256_t::operator|(const uint256_t &other) const {
    return uint256_t(high | other.high, low | other.low);
}

uint256_t& uint256_t::operator|=(const uint256_t &other) {
    high |= other.high;
    low |= other.low;
    return *this;
}

uint256_t uint256_t::operator^(const uint256_t &other) const {
    return uint256_t(high ^ other.high, low ^ other.low);
}

uint256_t& uint256_t::operator^=(const uint256_t &other) {
    high ^= other.high;
    low ^= other.low;
    return *this;
}

uint256_t uint256_t::operator~() const {
    return uint256_t(~high, ~low);
}

uint256_t uint256_t::operator-() const {
    uint256_t result;
    result.high = ~high;
    result.low = ~low;
    result.low += 1;
    
    if (result.low == 0)
        result.high += 1;

    return result;
}

bool uint256_t::operator==(const uint256_t &other) const {
    return high == other.high && low == other.low;
}

bool uint256_t::operator!=(const uint256_t &other) const {
    return !(*this == other);
}

bool uint256_t::operator!() const {
    return high == 0 && low == 0;
}

int uint256_t::popcount() const {
    int count = 0;
    count += __builtin_popcountll(high) + __builtin_popcountll(high >> 64);
    count += __builtin_popcountll(low) + __builtin_popcountll(low >> 64);
    return count;
}

uint256_t::operator bool() const {
    return high != 0 || low != 0;
}

uint256_t::operator int() const {
    return int(low);
}

uint256_t::operator uint8_t() const {
    return uint8_t(low);
}

uint256_t::operator uint64_t() const {
    return uint64_t(low);
}
