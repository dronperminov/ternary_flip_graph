#include "utils.h"

std::string prettyInt(size_t value) {
    std::stringstream ss;

    if (value < 1000)
        ss << value;
    else if (value < 1000000)
        ss << std::setprecision(2) << std::fixed << (value / 1000.0) << "K";
    else
        ss << std::setprecision(1) << std::fixed << (value / 1000000.0) << "M";

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
