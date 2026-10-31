#include "invariants_builder.h"

InvariantsBuilder::InvariantsBuilder(const std::vector<Ranks> &ranks) {
    this->ranks = ranks;
}

std::string InvariantsBuilder::getType() const {
    std::vector<Ranks> uniqueRanks;
    std::unordered_map<Ranks, int> powers;

    for (const Ranks& triplet : ranks) {
        if (powers.find(triplet) == powers.end())
            uniqueRanks.push_back(triplet);

        powers[triplet]++;
    }

    std::sort(uniqueRanks.begin(), uniqueRanks.end());
    std::stringstream type;
    bool printed = false;

    for (const Ranks& triplet : uniqueRanks) {
        if (printed)
            type << "+";

        int power = powers[triplet];
        if (power > 1)
            type << power;

        type << triplet;
        printed = true;
    }

    return type.str();
}

std::string InvariantsBuilder::getMinType() const {
    std::vector<std::vector<std::tuple<int, int, int>>> patterns(6);

    for (const Ranks& triplet : ranks) {
        int u = triplet.u;
        int v = triplet.v;
        int w = triplet.w;

        patterns[0].push_back({u, v, w});
        patterns[1].push_back({u, w, v});
        patterns[2].push_back({v, u, w});
        patterns[3].push_back({v, w, u});
        patterns[4].push_back({w, u, v});
        patterns[5].push_back({w, v, u});
    }

    size_t imin = 0;

    for (size_t i = 0; i < 6; i++)
        std::sort(patterns[i].begin(), patterns[i].end());

    for (size_t i = 1; i < 6; i++)
        if (patterns[i] < patterns[imin])
            imin = i;

    std::vector<Ranks> uniqueRanks;
    std::unordered_map<Ranks, int> powers;

    for (const auto& pattern : patterns[imin]) {
        Ranks triplet = {std::get<0>(pattern), std::get<1>(pattern), std::get<2>(pattern)};

        if (powers.find(triplet) == powers.end())
            uniqueRanks.push_back(triplet);

        powers[triplet]++;
    }

    std::sort(uniqueRanks.begin(), uniqueRanks.end());
    std::reverse(uniqueRanks.begin(), uniqueRanks.end());
    std::stringstream type;
    bool printed = false;

    for (const Ranks& triplet : uniqueRanks) {
        if (printed)
            type << "+";

        int power = powers[triplet];
        if (power > 1)
            type << power;

        type << triplet;
        printed = true;
    }

    return type.str();
}

std::string InvariantsBuilder::getF() const {
    std::unordered_map<Ranks, int> powers;

    for (const Ranks& triplet : ranks) {
        int u = triplet.u;
        int v = triplet.v;
        int w = triplet.w;

        powers[{u, v, w}]++;
        powers[{u, w, v}]++;
        powers[{v, u, w}]++;
        powers[{v, w, u}]++;
        powers[{w, u, v}]++;
        powers[{w, v, u}]++;
    }

    std::vector<Ranks> triplets;
    for (const auto& pair : powers)
        triplets.push_back(pair.first);

    std::sort(triplets.begin(), triplets.end(), [&](const Ranks& r1, const Ranks& r2) {
        int uvw1 = r1.u + r1.v + r1.w;
        int uvw2 = r2.u + r2.v + r2.w;
        if (uvw1 != uvw2)
            return uvw1 > uvw2;

        int uv1 = r1.u + r1.v;
        int uv2 = r2.u + r2.v;
        if (uv1 != uv2)
            return uv1 > uv2;

        return r1.u > r2.u;
    });

    std::stringstream type;
    bool printed = false;

    for (const Ranks& triplet : triplets) {
        if (printed)
            type << "+";

        int power = powers[triplet];
        if (power > 1)
            type << power;

        type << triplet;
        printed = true;
    }

    return type.str();
}
