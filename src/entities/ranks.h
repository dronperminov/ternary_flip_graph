#pragma once

#include <iostream>

struct Ranks {
    int u;
    int v;
    int w;

    bool operator==(const Ranks &ranks) const;
    bool operator<(const Ranks &ranks) const;

    friend std::ostream& operator<<(std::ostream &os, const Ranks &ranks);
};

namespace std {
    template <>
    struct hash<Ranks> {
        std::size_t operator()(const Ranks& ranks) const noexcept {
            std::size_t hu = std::hash<int>{}(ranks.u);
            std::size_t hv = std::hash<int>{}(ranks.v);
            std::size_t hw = std::hash<int>{}(ranks.w);

            return hu ^ hv ^ hw; 
        }
    };
}
