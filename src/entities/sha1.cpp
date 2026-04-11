#include "sha1.h"

std::string SHA1::get(const std::string &s) const {
    uint32_t digest[5] = {0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0};
    uint64_t transforms = 0;

    std::string buffer = "";
    std::istringstream is(s);

    while (1) {
        char sbuf[SHA1_BLOCK_BYTES];

        is.read(sbuf, SHA1_BLOCK_BYTES);
        buffer.append(sbuf, (size_t)is.gcount());

        if (buffer.size() != SHA1_BLOCK_BYTES)
            break;

        uint32_t block[SHA1_BLOCK_SIZE];
        splitToBlock(buffer, block);

        transform(digest, block);
        transforms++;

        buffer.clear();
    }

    return final(digest, buffer, transforms);
}

uint32_t SHA1::rol(uint32_t value, size_t bits) const {
    return (value << bits) | (value >> (32 - bits));
}

void SHA1::transform(uint32_t *digest, uint32_t *block) const {
    uint32_t a = digest[0];
    uint32_t b = digest[1];
    uint32_t c = digest[2];
    uint32_t d = digest[3];
    uint32_t e = digest[4];

    for (int i = 0; i < 80; i++) {
        if (i >= 16)
            block[i & 15] = rol(block[(i + 13) & 15] ^ block[(i + 8) & 15] ^ block[(i + 2) & 15] ^ block[i & 15], 1);

        if (i < 20) {
            e += ((b & (c ^ d)) ^ d) + 0x5A827999;
        }
        else if (i < 40) {
            e += (b ^ c ^ d) + 0x6ED9EBA1;
        }
        else if (i < 60) {
            e += (((b | c) & d) | (b & c)) + 0x8F1BBCDC;
        }
        else {
            e += (b ^ c ^ d) + 0xCA62C1D6;
        }

        uint32_t tmp = e + block[i & 15] + rol(a, 5);
        e = d;
        d = c;
        c = rol(b, 30);
        b = a;
        a = tmp;
    }

    digest[0] += a;
    digest[1] += b;
    digest[2] += c;
    digest[3] += d;
    digest[4] += e;
}

void SHA1::splitToBlock(const std::string &buffer, uint32_t *block) const {
    for (int i = 0; i < SHA1_BLOCK_SIZE; i++) {
        block[i] = 0;

        block[i] += (buffer[4 * i + 0] & 255) << 24;
        block[i] += (buffer[4 * i + 1] & 255) << 16;
        block[i] += (buffer[4 * i + 2] & 255) << 8;
        block[i] += (buffer[4 * i + 3] & 255) << 0;
    }
}

std::string SHA1::final(uint32_t *digest, std::string &buffer, uint64_t &transforms) const {
    uint64_t totalBits = (transforms * SHA1_BLOCK_BYTES + buffer.size()) * 8;

    buffer += (char) 0x80;

    size_t size = buffer.size();

    while (buffer.size() < SHA1_BLOCK_BYTES)
        buffer += (char) 0;

    uint32_t block[SHA1_BLOCK_SIZE];
    splitToBlock(buffer, block);

    if (size > SHA1_BLOCK_BYTES - 8) {
        transform(digest, block);
        transforms++;

        for (size_t i = 0; i < SHA1_BLOCK_SIZE - 2; i++)
            block[i] = 0;
    }

    block[SHA1_BLOCK_SIZE - 1] = (uint32_t) totalBits;
    block[SHA1_BLOCK_SIZE - 2] = (uint32_t) (totalBits >> 32);

    transform(digest, block);

    std::ostringstream result;

    for (size_t i = 0; i < 5; i++)
        result << std::hex << std::setfill('0') << std::setw(8) << digest[i];

    return result.str();
}
