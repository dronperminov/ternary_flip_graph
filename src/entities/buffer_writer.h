#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class BufferWriter {
    std::string path;
    std::vector<std::string> buffer;
    size_t bufferSize;
public:
    BufferWriter(const std::string &path, size_t bufferSize);

    void add(const std::string &line);

    ~BufferWriter();
private:
    void dump();
};
