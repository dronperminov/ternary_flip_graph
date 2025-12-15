#include "flip_set.h"

size_t FlipSet::size() const {
    return pairs.size();
}

void FlipSet::add(uint32_t index1, uint32_t index2) {
    pairs.push_back((index1 << 16) | index2);
}

void FlipSet::remove(uint32_t index1, uint32_t index2) {
    uint32_t pair1 = (index1 << 16) | index2;
    uint32_t pair2 = (index2 << 16) | index1;

    for (size_t i = 0; i < pairs.size(); i++) {
        if (pairs[i] == pair1 || pairs[i] == pair2) {
            pairs[i] = pairs.back();
            pairs.pop_back();
            return;
        }
    }
}

void FlipSet::remove(uint32_t index) {
    for (size_t i = 0; i < pairs.size(); i++) {
        if (index1(i) == index || index2(i) == index) {
            pairs[i--] = pairs.back();
            pairs.pop_back();
        }
    }
}

bool FlipSet::contains(uint32_t index1, uint32_t index2) const {
    uint32_t pair1 = (index1 << 16) | index2;
    uint32_t pair2 = (index2 << 16) | index1;

    for (size_t i = 0; i < pairs.size(); i++)
        if (pairs[i] == pair1 || pairs[i] == pair2)
            return true;

    return false;
}

void FlipSet::clear() {
    pairs.clear();
}

uint32_t FlipSet::index1(size_t i) const {
    return pairs[i] >> 16;
}

uint32_t FlipSet::index2(size_t i) const {
    return pairs[i] & 0xFFFF;
}
