#include "invariants_builder.h"

InvariantsBuilder::InvariantsBuilder(const std::vector<Ranks> &ranks) {
    this->ranks = ranks;
}

std::string InvariantsBuilder::getType() const {
    std::unordered_map<Ranks, int> powers;
    for (const Ranks& rankTriplet : ranks)
        powers[rankTriplet]++;

    std::vector<std::string> types;
    for (const auto &power: powers) {
        std::stringstream type;

        if (power.second > 1)
            type << power.second;

        type << power.first;
        types.push_back(type.str());
    }

    std::sort(types.begin(), types.end());
    std::stringstream type;

    for (size_t i = 0; i < types.size(); i++) {
        if (i > 0)
            type << "+";

        type << types[i];
    }

    return type.str();
}
