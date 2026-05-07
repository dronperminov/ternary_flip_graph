#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <filesystem>
#include "algebra/fraction.h"

std::string prettyInt(size_t value);
std::string prettyTime(double elapsed);

size_t parseNatural(std::string value);
Fraction parseFraction(const std::string &value);
std::vector<Fraction> parseFractions(const std::string &values);

bool makeDirectory(const std::string &path);
bool endsWith(const std::string &s, const std::string &substr);
int getMaxMatrixElements(const std::string &path, bool multiple);
std::vector<std::mt19937> initRandomGenerators(int seed, int count);
std::string getDimension(int n1, int n2, int n3, bool sorted = false);
