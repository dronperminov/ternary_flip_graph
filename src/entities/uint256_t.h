#pragma once

#include <iostream>
#include <cinttypes>

class uint256_t {
    __uint128_t high;
    __uint128_t low;
public:
    uint256_t();
    uint256_t(int value);
    uint256_t(unsigned int value);
    uint256_t(__uint128_t high, __uint128_t low);

    uint256_t operator<<(int n) const;
    uint256_t operator>>(int n) const;

    uint256_t operator&(const uint256_t &other) const;
    uint256_t& operator&=(const uint256_t &other);

    uint256_t operator&(int value) const;

    uint256_t operator|(const uint256_t &other) const;
    uint256_t& operator|=(const uint256_t &other);

    uint256_t operator^(const uint256_t &other) const;
    uint256_t& operator^=(const uint256_t &other);

    uint256_t operator~() const;
    uint256_t operator-() const;

    bool operator==(const uint256_t &other) const;
    bool operator!=(const uint256_t &other) const;

    bool operator!() const;

    int popcount() const;

    operator bool() const;
    operator int() const;
    operator uint8_t() const;
    operator uint64_t() const;
};
