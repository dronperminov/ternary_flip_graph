#pragma once

#include <unordered_map>
#include <string>

extern const std::unordered_map<std::string, int> KNOWN_RANKS_Q;
extern const std::unordered_map<std::string, int> KNOWN_RANKS_Z;
extern const std::unordered_map<std::string, int> KNOWN_RANKS_ZT;
extern const std::unordered_map<std::string, std::unordered_map<std::string, int>> KNOWN_RANKS;