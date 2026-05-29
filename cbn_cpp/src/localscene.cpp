#include "cbnetwork/localscene.hpp"
#include <iostream>

namespace cbnetwork {

void LocalAttractor::show() const {
    std::cout << "Network Index: " << network_index
              << ", Input Signal Index: [";
    for (size_t i = 0; i < relation_index.size(); ++i) {
        std::cout << relation_index[i] << (i == relation_index.size() - 1 ? "" : ",");
    }
    std::cout << "], Scene: " << local_scene
              << ", Global Index: " << g_index
              << ", Local Index: " << l_index
              << ", States: ";
    for (const auto& state : l_states) {
        std::cout << "[";
        for (size_t i = 0; i < state->l_variable_values.size(); ++i) {
            std::cout << state->l_variable_values[i] << (i == state->l_variable_values.size() - 1 ? "" : ",");
        }
        std::cout << "]";
    }
    std::cout << std::endl;
}

void LocalAttractor::show_short() const {
    std::cout << "Net. Index: " << network_index
              << ", Attrac. Index: " << l_index
              << ", States: ";
    for (const auto& state : l_states) {
        std::cout << "[";
        for (size_t i = 0; i < state->l_variable_values.size(); ++i) {
            std::cout << state->l_variable_values[i] << (i == state->l_variable_values.size() - 1 ? "" : ",");
        }
        std::cout << "]";
    }
    std::cout << std::endl;
}

} // namespace cbnetwork
