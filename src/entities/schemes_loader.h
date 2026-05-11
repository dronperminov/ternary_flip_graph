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
    std::vector<std::string> extensions;
public:
    SchemesLoader(const std::vector<std::string> &directories);

    std::vector<FractionalScheme> load(int n1, int n2, int n3, int rank, bool verify);
private:
    void loadFromDirectory(const std::string &directory, std::vector<FractionalScheme> &schemes, int n1, int n2, int n3, int rank, bool verify);
};
