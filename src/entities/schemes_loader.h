#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <filesystem>

#include "../utils.h"
#include "../schemes/fractional_scheme.h"

class SchemesLoader {
    std::unordered_map<std::string, std::vector<FractionalScheme>> dimension2schemes;
    std::vector<std::string> directories;
    std::unordered_map<std::string, std::vector<std::string>> ring2extensions;
public:
    SchemesLoader(const std::vector<std::string> &directories);

    std::vector<FractionalScheme> load(int n1, int n2, int n3, int rank, const std::string& ring, size_t maxCount, bool verify);
private:
    void loadFromDirectory(const std::string &directory, std::vector<FractionalScheme> &schemes, int n1, int n2, int n3, int rank, const std::string& ring, size_t maxCount, bool verify);
};
