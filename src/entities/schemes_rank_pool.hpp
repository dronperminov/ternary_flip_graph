#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include "../utils.h"
#include "schemes_pool.hpp"

template <typename Scheme>
class SchemesRankPool {
    std::string dimension;
    size_t maxSize;
    bool uniqueOnly;
    std::string path;
    std::string format;

    std::vector<int> ranks;
    std::unordered_map<int, SchemesPool<Scheme>> rank2pool;
public:
    SchemesRankPool(const std::string &dimension, size_t maxSize, bool uniqueOnly, const std::string &path, const std::string &format);

    int minRank() const;
    int maxRank() const;

    size_t minRankSize() const;
    size_t size(int rank) const;

    double minFillRatio() const;
    double fillRatio(int rank) const;

    bool add(const Scheme &scheme, bool save = false);
    void copyRandom(Scheme &scheme, std::mt19937 &generator) const;
    void copyRandomMinRank(Scheme &scheme, std::mt19937 &generator) const;
    void print(int knownRank) const;
private:
    int getRandomRank(std::mt19937 &generator, double alpha) const;
};

template <typename Scheme>
SchemesRankPool<Scheme>::SchemesRankPool(const std::string &dimension, size_t maxSize, bool uniqueOnly, const std::string &path, const std::string &format) {
    this->dimension = dimension;
    this->maxSize = maxSize;
    this->uniqueOnly = uniqueOnly;
    this->path = path;
    this->format = format;
}

template <typename Scheme>
int SchemesRankPool<Scheme>::minRank() const {
    if (ranks.empty())
        return -1;

    return ranks[0];
}

template <typename Scheme>
int SchemesRankPool<Scheme>::maxRank() const {
    if (ranks.empty())
        return -1;

    return ranks.back();
}

template <typename Scheme>
size_t SchemesRankPool<Scheme>::minRankSize() const {
    if (ranks.empty())
        return 0;

    return rank2pool.at(ranks[0]).size();
}

template <typename Scheme>
size_t SchemesRankPool<Scheme>::size(int rank) const {
    if (rank2pool.find(rank) == rank2pool.end())
        return 0;

    return rank2pool.at(rank).size();
}

template <typename Scheme>
double SchemesRankPool<Scheme>::minFillRatio() const {
    if (ranks.empty())
        return 0;

    return rank2pool.at(ranks[0]).size() / double(maxSize);
}

template <typename Scheme>
double SchemesRankPool<Scheme>::fillRatio(int rank) const {
    return size(rank) / double(maxSize);
}

template <typename Scheme>
bool SchemesRankPool<Scheme>::add(const Scheme &scheme, bool save) {
    int rank = scheme.getRank();

    if (rank2pool.find(rank) == rank2pool.end()) {
        std::stringstream ss;
        ss << path << "/rank" << rank;
        std::string rankPath = ss.str();

        rank2pool.emplace(rank, SchemesPool<Scheme>(maxSize, uniqueOnly, rankPath, format));
        ranks.push_back(rank);
        std::sort(ranks.begin(), ranks.end());
    }

    return rank2pool.at(rank).add(scheme, save);
}

template <typename Scheme>
void SchemesRankPool<Scheme>::print(int knownRank) const {
    std::cout << "+-----------+------+---------+" << std::endl;

    for (size_t i = 0; i < ranks.size(); i++) {
        std::cout << "| ";
        std::cout << std::setw(9) << (i == 0 ? dimension : "") << " | ";
        std::cout << std::setw(4) << ranks[i] << " | ";
        std::cout << std::setw(7) << rank2pool.at(ranks[i]).size();
        std::cout << " |";

        if (i == 0) {
            if (ranks[0] < knownRank)
                std::cout << " improved (known: " << knownRank << ")";
            else if (ranks[0] > knownRank)
                std::cout << " known: " << knownRank;
        }

        std::cout << std::endl;
    }
}

template <typename Scheme>
int SchemesRankPool<Scheme>::getRandomRank(std::mt19937 &generator, double alpha) const {
    std::vector<double> weights(ranks.size());
    weights[0] = 1.0;

    double total = 1.0;

    for (size_t i = 1; i < ranks.size(); i++) {
        weights[i] = weights[i - 1] * alpha;
        total += weights[i];
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double randomWeight = uniform(generator) * total;
    double sum = 0;

    for (size_t i = 0; i < ranks.size(); i++) {
        sum += weights[i];

        if (randomWeight <= sum)
            return ranks[i];
    }

    return ranks.back();
}

template <typename Scheme>
void SchemesRankPool<Scheme>::copyRandom(Scheme &scheme, std::mt19937 &generator) const {
    int rank = getRandomRank(generator, 0.7);
    rank2pool.at(rank).copyRandom(scheme, generator);
}

template <typename Scheme>
void SchemesRankPool<Scheme>::copyRandomMinRank(Scheme &scheme, std::mt19937 &generator) const {
    rank2pool.at(ranks[0]).copyRandom(scheme, generator);
}
