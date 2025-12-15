#pragma once

#include <iostream>
#include <vector>
#include <cstdint>

class FlipSet {
    std::vector<uint32_t> pairs;
public:
    size_t size() const;

    void add(uint32_t index1, uint32_t index2);
    void remove(uint32_t index1, uint32_t index2);
    void remove(uint32_t index);
    bool contains(uint32_t index1, uint32_t index2) const;
    void clear();

    uint32_t index1(size_t i) const;
    uint32_t index2(size_t i) const;
};
