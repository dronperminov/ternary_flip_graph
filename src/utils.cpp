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