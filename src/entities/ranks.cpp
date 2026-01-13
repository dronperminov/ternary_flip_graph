#include "ranks.h"

bool Ranks::operator==(const Ranks &ranks) const {
    return u == ranks.u && v == ranks.v && w == ranks.w;
}

bool Ranks::operator<(const Ranks &ranks) const {
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
