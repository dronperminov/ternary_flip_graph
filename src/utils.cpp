#include "utils.h"

std::string prettyInt(size_t value) {
    std::stringstream ss;

    if (value < 1000)
        ss << value;
    else if (value < 1000000)
        ss << std::setprecision(2) << std::fixed << (value / 1000.0) << "K";
    else if (value < 1000000000)
        ss << std::setprecision(1) << std::fixed << (value / 1000000.0) << "M";
    else
        ss << std::setprecision(1) << std::fixed << (value / 1000000000.0) << "B";

    return ss.str();
}

std::string prettyTime(double elapsed) {
    std::stringstream ss;

    if (elapsed < 60) {
        ss << std::setprecision(2) << std::fixed << elapsed;
    }
    else {
        int seconds = int(elapsed + 0.5);
        int hours = seconds / 3600;
        int minutes = (seconds % 3600) / 60;

        ss << std::setw(2) << std::setfill('0') << hours << ":";
        ss << std::setw(2) << std::setfill('0') << minutes << ":";
        ss << std::setw(2) << std::setfill('0') << (seconds % 60);
    }

    return ss.str();
}

size_t parseNatural(std::string value) {
    size_t multiplier = 1;

    if (value.back() == 'K' || value.back() == 'k') {
        multiplier = 1000;
    }
    else if (value.back() == 'M' || value.back() == 'm') {
        multiplier = 1000000;
    }
    else if (value.back() == 'B' || value.back() == 'b') {
        multiplier = 1000000000;
    }

    if (multiplier > 1)
        value.pop_back();

    return std::stoul(value) * multiplier;
}

bool makeDirectory(const std::string &path) {
    std::error_code err;
    if (std::filesystem::create_directories(path, err))
        return true;

    if (std::filesystem::exists(path))
        return true;

    std::cerr << "Unable to create directory \"" << path << "\": " << err.message() << std::endl;
    return false;
}

int getMaxMatrixElements(const std::string &path, bool multiple) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Unable to open file \"" << path << "\"" << std::endl;
        return -1;
    }

    int count = 1;
    if (multiple)
        f >> count;

    int n1, n2, n3;
    f >> n1 >> n2 >> n3;
    f.close();

    int maxMatrixElements = std::max(n1 * n2, std::max(n2 * n3, n3 * n1));

    if (maxMatrixElements > 128) {
        std::cerr << "Max matrix elements too big (> 128): " << n1 << "x" << n2 << "x" << n3 << std::endl;
        return -1;
    }

    return maxMatrixElements;
}
