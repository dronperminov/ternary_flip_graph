#include "buffer_writer.h"

BufferWriter::BufferWriter(const std::string &path, size_t bufferSize) {
    this->path = path;
    this->bufferSize = bufferSize;
    buffer.reserve(bufferSize);
}

void BufferWriter::add(const std::string &line) {
    buffer.push_back(line);

    if (buffer.size() >= bufferSize)
        dump();
}

BufferWriter::~BufferWriter() {
    dump();
}

void BufferWriter::dump() {
    if (buffer.empty())
        return;

    std::ofstream f(path, std::ios::app);

    for (const std::string &line : buffer)
        f << line << std::endl;

    f.close();
    buffer.clear();
}
