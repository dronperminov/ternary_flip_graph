#include "schemes_loader.h"

SchemesLoader::SchemesLoader(const std::vector<std::string> &directories) {
    this->directories = directories;
    this->ring2extensions = {
        {"ZT", {"ZT.txt"}},
        {"Z", {"ZT.txt", "Z.txt"}},
        {"Q", {"ZT.txt", "Z.txt", "Q.txt"}}
    };
}

std::vector<FractionalScheme> SchemesLoader::load(int n1, int n2, int n3, int rank, const std::string& ring, size_t maxCount, bool verify) {
    std::string dimension = getDimension(n1, n2, n3);

    if (std::min(n1, std::min(n2, n3)) == 1)
        return {FractionalScheme(n1, n2, n3)};

    if (dimension2schemes.find(dimension) != dimension2schemes.end())
        return dimension2schemes.at(dimension);

    std::string key = getDimension(n1, n2, n3, true);
    std::vector<FractionalScheme> schemes;
    for (const std::string &directory : directories)
        loadFromDirectory(directory + "/" + key + "/rank" + std::to_string(rank), schemes, n1, n2, n3, rank, ring, maxCount, verify);

    dimension2schemes[dimension] = schemes;
    std::cout << "Loader read " << schemes.size() << " schemes (" << dimension << ": " << rank << ")" << std::endl;
    return schemes;
}

void SchemesLoader::loadFromDirectory(const std::string &directory, std::vector<FractionalScheme> &schemes, int n1, int n2, int n3, int rank, const std::string& ring, size_t maxCount, bool verify) {
    if (!std::filesystem::exists(directory))
        return;

    std::vector<std::string> extensions = ring2extensions.at(ring);

    for (auto it = std::filesystem::directory_iterator(directory); it != std::filesystem::directory_iterator(); it++) {
        if (!it->is_regular_file())
            continue;

        std::string path = it->path().string();
        if (!endsWith(path, extensions))
            continue;

        FractionalScheme scheme;
        if (scheme.read(path, verify, !endsWith(path, "Q.txt")) && scheme.getRank() == rank && scheme.setDimension(n1, n2, n3))
            schemes.emplace_back(scheme);

        if (maxCount && schemes.size() >= maxCount)
            break;
    }
}
