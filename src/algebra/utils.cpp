#include "utils.h"

int64_t modInverse(int64_t value, int64_t mod) {
    int64_t x0 = 1, x1 = 0;
    int64_t r0 = value, r1 = mod;

    while (r1 != 0) {
        int64_t q = r0 / r1;
        int64_t r2 = r0 - q * r1;
        int64_t x2 = x0 - q * x1;

        r0 = r1;
        r1 = r2;
        x0 = x1;
        x1 = x2;
    }

    if (r0 != 1)
        throw std::runtime_error("mod inverse does not exists");

    if (x0 < 0)
        x0 += mod;

    return x0;
}
