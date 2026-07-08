#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <vector>
#include <sstream>
#include <string>
#include <chrono>
#include <algorithm>
#include <filesystem>
#include "algebra/fraction.h"

std::string prettyInt(size_t value);
std::string prettyTime(double elapsed);
std::string prettyTime(const std::chrono::high_resolution_clock::time_point& t1, const std::chrono::high_resolution_clock::time_point& t2);

size_t parseNatural(std::string value);
Fraction parseFraction(const std::string &value);
std::vector<Fraction> parseFractions(const std::string &values);

bool makeDirectory(const std::string &path);
bool endsWith(const std::string &s, const std::string &substr);
bool endsWith(const std::string &s, const std::vector<std::string> &substrs);
std::string join(const std::vector<std::string> &values, const std::string& delimeter = "");
int getMaxMatrixElements(const std::string &path, bool multiple);
std::vector<std::mt19937> initRandomGenerators(int seed, int count);

std::string getDimension(int n1, int n2, int n3, bool sorted = false);
int digitsCount(size_t n);
bool isPowerOfTwo(int n);

std::vector<std::string> getSchemePathsFromDirectory(const std::string &inputPath, const std::vector<std::string> &extensions);
std::vector<std::string> getSchemePathsFromDirectoryRecursive(const std::string &inputPath, const std::vector<std::string> &extensions);
std::vector<std::string> getSchemePathsFromFile(const std::string &inputPath, const std::vector<std::string> &extensions);
std::vector<std::string> getSchemePathsFromPath(const std::string &inputPath, const std::vector<std::string> &extensions, bool recursive);
