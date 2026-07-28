#include "invariants_builder.h"

InvariantsBuilder::InvariantsBuilder(const std::vector<Ranks> &ranks) {
    this->ranks = ranks;
}

std::string InvariantsBuilder::getType() const {
    std::vector<Ranks> uniqueRanks;
    std::unordered_map<Ranks, int> powers;

    for (const Ranks& rankTriplet : ranks) {
        if (powers.find(rankTriplet) == powers.end())
            uniqueRanks.push_back(rankTriplet);

        powers[rankTriplet]++;
    }

    std::sort(uniqueRanks.begin(), uniqueRanks.end());
    std::stringstream type;
    bool printed = false;

    for (const Ranks& rankTriplet : uniqueRanks) {
        if (printed)
            type << "+";

        int power = powers[rankTriplet];
        if (power > 1)
            type << power;

        type << rankTriplet;
        printed = true;
    }

    return type.str();
}
