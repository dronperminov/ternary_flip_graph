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
    bool uniqueOnly;
    std::string path;
    std::string format;
    size_t maxSave;
    size_t saved;
    std::vector<Scheme> schemes;
    size_t index;
    bool hasDirectory;

    int minComplexity;
    int maxComplexity;
    int minFlips;
    int maxFlips;

    std::unordered_set<std::string> hashes;
    SHA1 sha1;
public:
    SchemesPool(size_t maxSize, bool uniqueOnly, const std::string &path, const std::string &format, size_t maxSave = 32768);

    size_t size() const;
    int getMinComplexity() const;
    int getMaxComplexity() const;
    int getMinFlips() const;
    int getMaxFlips() const;

    bool add(const Scheme &scheme, bool save);
    void copyRandom(Scheme &scheme, std::mt19937 &generator) const;
private:
    void saveScheme(const Scheme &scheme);
};

template <typename Scheme>
SchemesPool<Scheme>::SchemesPool(size_t maxSize, bool uniqueOnly, const std::string &path, const std::string &format, size_t maxSave) {
    this->maxSize = maxSize;
    this->uniqueOnly = uniqueOnly;
    this->path = path;
    this->format = format;
    this->maxSave = maxSave;
    this->saved = 0;
    this->index = 0;
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
    if (uniqueOnly) {
        std::string hash = sha1.get(scheme.getHash());
        if (hashes.find(hash) != hashes.end())
            return false;

        hashes.insert(hash);
    }

    if (schemes.size() < maxSize)
        schemes.emplace_back(Scheme());

    schemes[index].copy(scheme);
    index = (index + 1) % maxSize;

    if (save && saved < maxSave)
        saveScheme(scheme);

    int complexity = scheme.getComplexity();
    int flips = scheme.getAvailableFlips();

    minComplexity = std::min(minComplexity, complexity);
    maxComplexity = std::max(maxComplexity, complexity);
    minFlips = std::min(minFlips, flips);
    maxFlips = std::max(maxFlips, flips);

    return true;
}

template <typename Scheme>
void SchemesPool<Scheme>::copyRandom(Scheme &scheme, std::mt19937 &generator) const {
    scheme.copy(schemes[generator() % schemes.size()]);
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

    saved++;
}
