#pragma once

#include <algorithm>
#include <sstream>
#include <vector>
#include <unordered_map>
#include "ranks.h"

class InvariantsBuilder {
    std::vector<Ranks> ranks;
public:
    InvariantsBuilder(const std::vector<Ranks> &ranks);

    std::string getType() const;
};
