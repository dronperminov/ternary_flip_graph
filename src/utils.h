#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

std::string prettyInt(size_t value);
std::string prettyTime(double elapsed);
size_t parseNatural(std::string value);