#ifndef PATH_UTILS_HPP
#define PATH_UTILS_HPP

#include <string>
#include <filesystem>

namespace cbnetwork {

namespace fs = std::filesystem;

/**
 * Finds the project root by looking for the .git directory or cbn_cpp folder.
 */
inline fs::path find_project_root() {
    fs::path curr = fs::current_path();
    while (curr.has_parent_path()) {
        if (fs::exists(curr / ".git") || fs::exists(curr / "cbn_cpp")) {
            return curr;
        }
        curr = curr.parent_path();
    }
    return fs::current_path(); // Fallback to current path
}

/**
 * Creates the dynamic output directory based on configuration.
 */
inline std::string create_output_directory(int topology, int n_networks, int n_vars, int n_samples) {
    fs::path root = find_project_root();
    fs::path base_dir = root / "phd_experiments";

    std::string folder_name = "output_" + std::to_string(topology) + "_" +
                              std::to_string(n_networks) + "_" +
                              std::to_string(n_vars) + "_" +
                              std::to_string(n_samples);

    fs::path output_path = base_dir / folder_name;

    if (!fs::exists(output_path)) {
        fs::create_directories(output_path);
    }

    return output_path.string();
}

} // namespace cbnetwork

#endif // PATH_UTILS_HPP
