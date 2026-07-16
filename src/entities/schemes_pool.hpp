#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <climits>
#include <unordered_set>

#include "../entities/sha1.h"

template <typename Scheme>
class SchemesPool {
    size_t maxSize;
    std::string uniqueType;
    std::string path;
    std::string format;

    std::vector<Scheme> schemes;
    std::vector<int> schemeFlips;
    size_t totalFlips;
    size_t changes;
    size_t index;
    bool hasDirectory;

    int minComplexity;
    int maxComplexity;
    int minFlips;
    int maxFlips;

    std::unordered_set<std::string> hashes;
    SHA1 sha1;
public:
    SchemesPool(size_t maxSize, const std::string uniqueType, const std::string &path, const std::string &format);

    size_t size() const;
    size_t getDiff() const;
    int getMinComplexity() const;
    int getMaxComplexity() const;
    int getMinFlips() const;
    int getMaxFlips() const;

    bool add(const Scheme &scheme, bool save);
    bool contains(const Scheme &scheme) const;
    void copyRandom(Scheme &scheme, std::mt19937 &generator, bool selectUniform = true) const;
    void resetDiff();
private:
    void saveScheme(const Scheme &scheme);
    std::string getHash(const Scheme &scheme) const;
};

template <typename Scheme>
SchemesPool<Scheme>::SchemesPool(size_t maxSize, const std::string uniqueType, const std::string &path, const std::string &format) {
    this->maxSize = maxSize;
    this->uniqueType = uniqueType;
    this->path = path;
    this->format = format;
    this->index = 0;
    this->totalFlips = 0;
    this->changes = 0;
    this->hasDirectory = false;

    this->minComplexity = INT_MAX;
    this->maxComplexity = 0;
    this->minFlips = INT_MAX;
    this->maxFlips = 0;
}

template <typename Scheme>
size_t SchemesPool<Scheme>::size() const {
    return schemes.size();
}

template <typename Scheme>
size_t SchemesPool<Scheme>::getDiff() const {
    return changes;
}

template <typename Scheme>
int SchemesPool<Scheme>::getMinComplexity() const {
    return minComplexity;
}

template <typename Scheme>
int SchemesPool<Scheme>::getMaxComplexity() const {
    return maxComplexity;
}

template <typename Scheme>
int SchemesPool<Scheme>::getMinFlips() const {
    return minFlips;
}

template <typename Scheme>
int SchemesPool<Scheme>::getMaxFlips() const {
    return maxFlips;
}

template <typename Scheme>
bool SchemesPool<Scheme>::add(const Scheme &scheme, bool save) {
    if (!uniqueType.empty()) {
        std::string hash = sha1.get(getHash(scheme));
        if (hashes.find(hash) != hashes.end())
            return false;

        hashes.insert(hash);
    }

    int complexity = scheme.getComplexity();
    int flips = scheme.getAvailableFlips();

    if (schemes.size() < maxSize) {
        schemes.emplace_back(Scheme());
        schemeFlips.push_back(flips);

        if (save)
            saveScheme(scheme);
    }
    else {
        totalFlips -= schemeFlips[index];
        schemeFlips[index] = flips;
    }

    schemes[index].copy(scheme);
    index = (index + 1) % maxSize;

    minComplexity = std::min(minComplexity, complexity);
    maxComplexity = std::max(maxComplexity, complexity);
    minFlips = std::min(minFlips, flips);
    maxFlips = std::max(maxFlips, flips);
    totalFlips += flips;
    changes++;

    return true;
}

template <typename Scheme>
bool SchemesPool<Scheme>::contains(const Scheme &scheme) const {
    std::string hash = getHash(scheme);
    return hashes.find(hash) != hashes.end();
}

template <typename Scheme>
void SchemesPool<Scheme>::copyRandom(Scheme &scheme, std::mt19937 &generator, bool selectUniform) const {
    if (totalFlips == 0 || selectUniform) {
        scheme.copy(schemes[generator() % schemes.size()]);
        return;
    }

    std::uniform_int_distribution<size_t> uniform(0, totalFlips);

    size_t randomWeight = uniform(generator);
    size_t sum = 0;

    for (size_t i = 0; i < schemes.size(); i++) {
        sum += schemeFlips[i];

        if (randomWeight <= sum) {
            scheme.copy(schemes[i]);
            return;
        }
    }

    scheme.copy(schemes.back());
}

template <typename Scheme>
void SchemesPool<Scheme>::resetDiff() {
    changes = 0;
}

template <typename Scheme>
void SchemesPool<Scheme>::saveScheme(const Scheme &scheme) {
    if (!hasDirectory && !std::filesystem::exists(path)) {
        makeDirectory(path);
        hasDirectory = true;
    }

    std::stringstream ss;
    ss << path << "/";
    ss << scheme.getDimension();
    ss << "_m" << scheme.getRank();
    ss << "_" << sha1.get(scheme.getHash());
    ss << "_" << scheme.getRing();
    ss << "." << format;

    if (format == "json")
        scheme.saveJson(ss.str());
    else
        scheme.saveTxt(ss.str());
}

template <typename Scheme>
std::string SchemesPool<Scheme>::getHash(const Scheme &scheme) const {
    if (uniqueType == "structure")
        return scheme.getStructureHash();

    if (uniqueType == "type")
        return scheme.getTypeInvariant();

    return scheme.getHash();
}
