#include "arg_parser.h"


Argument::Argument(const std::string &longName, const std::string &shortName, ArgType type, const std::string &description, const std::string &section, const std::vector<std::string> &choices, const std::string &defaultValue, bool required) {
    this->longName = longName;
    this->shortName = shortName;
    this->type = type;
    this->description = description;
    this->section = section;
    this->choices = choices;
    this->defaultValue = defaultValue;
    this->value = defaultValue;
    this->required = required;
    this->isSet = false;
}

std::string Argument::getHelpName() const {
    std::stringstream ss;

    if (longName.empty()) {
        ss << shortName;
    }
    else if (shortName.empty()) {
        ss << longName;
    }
    else {
        ss << shortName + ", " + longName;
    }

    ss << " " << getTypeHint();
    return ss.str();
}

std::string Argument::getHelpDescription() const {
    std::stringstream ss;
    ss << description;

    if (type != ArgType::Flag && !defaultValue.empty())
        ss << " (default: " << defaultValue << ")";

    if (required)
        ss << " [REQUIRED]";

    return ss.str();
}

std::string Argument::getTypeHint() const {
    if (!choices.empty()) {
        std::stringstream ss;
        ss << "{";

        for (size_t i = 0; i < choices.size(); i++)
            ss << (i > 0 ? ", " : "") << choices[i];

        ss << "}";
        return ss.str();
    }

    if (type == ArgType::Natural)
        return "INT";

    if (type == ArgType::Real)
        return "REAL";

    if (type == ArgType::String)
        return "STR";

    if (type == ArgType::Path)
        return "PATH";

    return "";
}

size_t Argument::getWidth() const {
    size_t typeHintWidth = 1 + getTypeHint().length();

    if (longName.empty())
        return shortName.length() + typeHintWidth;

    if (shortName.empty())
        return longName.length() + typeHintWidth;

    return longName.length() + 2 + shortName.length() + typeHintWidth;
}

bool Argument::validate(const std::string &parsedName) const {
    if (type != ArgType::Flag && !isValidChoice(value)) {
        std::cerr << "Invalid value for argument \"" << parsedName << "\": " << value << " is not valid choice. Valid choices are: ";

        for (size_t i = 0; i < choices.size(); i++)
            std::cerr << (i > 0 ? ", " : "") << choices[i];

        std::cerr << std::endl;
        return false;
    }

    if (type == ArgType::String || type == ArgType::Path)
        return true;

    if (type == ArgType::Natural && !isNatural(value)) {
        std::cerr << "Invalid value for argument \"" << parsedName << "\": " << value << " is not natural" << std::endl;
        return false;
    }

    if (type == ArgType::Real && !isReal(value)) {
        std::cerr << "Invalid value for argument \"" << parsedName << "\": " << value << " is not real" << std::endl;
        return false;
    }

    return true;
}

bool Argument::isNatural(const std::string &value) const {
    size_t size = value.size();

    if (value.back() == 'K' || value.back() == 'k' || value.back() == 'M' || value.back() == 'm' || value.back() == 'B' || value.back() == 'b')
        size--;

    for (size_t i = 0; i < size; i++)
        if (value[i] < '0' || value[i] > '9')
            return false;

    return std::stoul(value) > 0;
}

bool Argument::isReal(const std::string &value) const {
    if (value.size() == 0)
        return false;

    size_t start = value[0] == '-' ? 1 : 0;
    bool point = false;

    for (size_t i = start; i < value.size(); i++) {
        if (value[i] == '.') {
            if (point)
                return false;

            point = true;
        }
        else if (value[i] < '0' || value[i] > '9')
            return false;
    }

    return true;
}

bool Argument::isValidChoice(const std::string &value) const {
    if (choices.empty())
        return true;

    return std::find(choices.begin(), choices.end(), value) != choices.end();
}

ArgParser::ArgParser(const std::string &name, const std::string &description) {
    this->name = name;
    this->description = description;
    this->maxArgWidth = 30;

    addSection("Main options");
    add("--help", "-h", ArgType::Flag, "Show this help message");
}

void ArgParser::addSection(const std::string &section) {
    sections.push_back(section);
}

void ArgParser::addChoices(const std::string &longName, const std::string &shortName, ArgType type, const std::string &description, const std::vector<std::string> &choices, const std::string &defaultValue, bool required) {
    Argument argument(longName, shortName, type, description, sections.back(), choices, defaultValue, required);
    arguments.emplace_back(argument);

    if (!longName.empty())
        arg2index[longName] = arguments.size() - 1;

    if (!shortName.empty())
        arg2index[shortName] = arguments.size() - 1;

    maxArgWidth = std::max(maxArgWidth, argument.getWidth());
}

void ArgParser::addChoices(const std::string &name, ArgType type, const std::string &description, const std::vector<std::string> &choices, const std::string &defaultValue, bool required) {
    if (name.find("--") == 0) {
        addChoices(name, "", type, description, choices, defaultValue, required);
    }
    else if (name.find("-") == 0) {
        addChoices("", name, type, description, choices, defaultValue, required);
    }
    else
        throw std::runtime_error("Argument name must starts with \"-\" or \"--\"");
}

void ArgParser::add(const std::string &longName, const std::string &shortName, ArgType type, const std::string &description, const std::string &defaultValue, bool required) {
    addChoices(longName, shortName, type, description, {}, defaultValue, required);
}

void ArgParser::add(const std::string &name, ArgType type, const std::string &description, const std::string &defaultValue, bool required) {
    addChoices(name, type, description, {}, defaultValue, required);
}

bool ArgParser::parse(int argc, char *argv[]) {
    if (argc > 0)
        name = argv[0];

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            help();
            return false;
        }

        auto it = arg2index.find(arg);
        if (it == arg2index.end()) {
            std::cerr << "Unknown argument \"" << arg << "\"" << std::endl;
            return false;
        }

        Argument& argument = arguments[it->second];
        if (argument.isSet) {
            std::cerr << "Argument \"" << arg << "\" has already been set" << std::endl;
            return false;
        }

        argument.isSet = true;

        if (argument.type == ArgType::Flag) {
            argument.value = "true";
            continue;
        }

        if (i + 1 >= argc) {
            std::cerr << "Missing value for argument \"" << arg << "\"" << std::endl;
            return false;
        }

        argument.value = argv[++i];
        if (!argument.validate(arg))
            return false;
    }

    for (const auto& argument : arguments) {
        if (argument.required && !argument.isSet) {
            std::cerr << "Required argument \"" << (argument.longName.empty() ? argument.shortName : argument.longName) << "\" is missing" << std::endl;
            return false;
        }
    }

    return true;
}

bool ArgParser::isSet(const std::string &name) const {
    auto it = arg2index.find(name);
    if (it == arg2index.end())
        return false;

    return arguments[it->second].isSet;
}

void ArgParser::help() const {
    std::cout << description << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: " << name << " [ARGS...]" << std::endl;

    std::unordered_map<std::string, std::vector<const Argument *>> section2args;

    for (const auto &argument : arguments)
        section2args[argument.section].push_back(&argument);

    for (const auto &section : sections) {
        if (section2args.at(section).empty())
            continue;

        std::cout << std::endl << section << ":" << std::endl;

        for (const auto &argument : section2args.at(section)) {
            std::cout << std::left << std::setw(maxArgWidth + 2) << argument->getHelpName();
            std::cout << argument->getHelpDescription();
            std::cout << std::endl;
        }
    }
}

std::string ArgParser::get(const std::string &name) const {
    auto it = arg2index.find(name);
    if (it == arg2index.end())
        throw std::runtime_error("Argument not found: " + name);

    return arguments[it->second].value;
}

std::string ArgParser::operator[](const std::string &name) const {
    return get(name);
}
