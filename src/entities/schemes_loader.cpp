#include "schemes_loader.h"

SchemesLoader::SchemesLoader(const std::vector<std::string> &directories) {
    this->directories = directories;
    this->extensions = {"ZT.txt", "Z.txt", "Q.txt"};
}

std::vector<FractionalScheme> SchemesLoader::load(int n1, int n2, int n3, int rank, bool verify) {
    std::string dimension = getDimension(n1, n2, n3);

    if (std::min(n1, std::min(n2, n3)) == 1)
        return {FractionalScheme(n1, n2, n3)};

    if (dimension2schemes.find(dimension) != dimension2schemes.end())
        return dimension2schemes.at(dimension);

    std::string key = getDimension(n1, n2, n3, true);
    std::vector<FractionalScheme> schemes;
    for (const std::string &directory : directories)
        loadFromDirectory(directory + "/" + key + "/rank" + std::to_string(rank), schemes, n1, n2, n3, rank, verify);

    dimension2schemes[dimension] = schemes;
    std::cout << "Loader read " << schemes.size() << " schemes (" << dimension << ": " << rank << ")" << std::endl;
    return schemes;
}

void SchemesLoader::loadFromDirectory(const std::string &directory, std::vector<FractionalScheme> &schemes, int n1, int n2, int n3, int rank, bool verify) {
    std::vector<std::string> paths;

    for (auto it = std::filesystem::directory_iterator(directory); it != std::filesystem::directory_iterator(); it++) {
        if (!it->is_regular_file())
            continue;

        std::string path = it->path().string();
        if (endsWith(path, extensions))
            paths.push_back(path);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < paths.size(); i++) {
        FractionalScheme scheme;
        if (scheme.read(paths[i], verify, !endsWith(paths[i], "Q.txt")) && scheme.getRank() == rank && scheme.setDimension(n1, n2, n3))
            schemes.emplace_back(scheme);
    }
}
