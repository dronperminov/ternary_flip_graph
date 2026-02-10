#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <filesystem>

std::string prettyInt(size_t value);
std::string prettyTime(double elapsed);
size_t parseNatural(std::string value);
bool makeDirectory(const std::string &path);
int getMaxMatrixElements(const std::string &path, bool multiple);