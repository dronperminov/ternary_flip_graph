#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdint>

const int SHA1_BLOCK_SIZE = 16;
const int SHA1_BLOCK_BYTES = 64;

class SHA1 {
public:
    std::string get(const std::string &s) const;
private:
    uint32_t rol(uint32_t value, size_t bits) const;
    void transform(uint32_t *digest, uint32_t *block) const;
    void splitToBlock(const std::string &buffer, uint32_t *block) const;
    std::string final(uint32_t *digest, std::string &buffer, uint64_t &transforms) const;
};
