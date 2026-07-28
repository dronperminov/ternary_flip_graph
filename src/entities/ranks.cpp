#include "ranks.h"

bool Ranks::operator==(const Ranks &ranks) const {
    return u == ranks.u && v == ranks.v && w == ranks.w;
}

bool Ranks::operator<(const Ranks &ranks) const {
    int uvw1 = u + v + w;
    int uvw2 = ranks.u + ranks.v + ranks.w;
    if (uvw1 != uvw2)
        return uvw1 < uvw2;

    int uv1 = u + v;
    int uv2 = ranks.u + ranks.v;
    if (uv1 != uv2)
        return uv1 < uv2;

    if (u != ranks.u)
        return u < ranks.u;

    if (v != ranks.v)
        return v < ranks.v;

    return w < ranks.w;
}

std::ostream& operator<<(std::ostream &os, const Ranks &ranks) {
    if (ranks.u > 0)
        os << "x";

    if (ranks.u > 1)
        os << "^" << ranks.u;

    if (ranks.v > 0)
        os << "y";

    if (ranks.v > 1)
        os << "^" << ranks.v;

    if (ranks.w > 0)
        os << "z";

    if (ranks.w > 1)
        os << "^" << ranks.w;

    return os;
}
