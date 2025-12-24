#include "base_scheme.h"

BaseScheme::BaseScheme() : boolDistribution(0, 1), ijkDistribution(0, 2) {

}

int BaseScheme::getRank() const {
    return rank;
}

int BaseScheme::getDimension(int index) const {
    return dimension[index];
}

std::string BaseScheme::getDimension() const {
    std::stringstream ss;
    ss << dimension[0] << "x" << dimension[1] << "x" << dimension[2];
    return ss.str();
}

int BaseScheme::getAvailableFlips() const {
    return flips[0].size() + flips[1].size() + flips[2].size();
}
