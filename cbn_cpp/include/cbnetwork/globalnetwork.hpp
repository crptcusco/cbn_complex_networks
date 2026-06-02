#ifndef GLOBALNETWORK_HPP
#define GLOBALNETWORK_HPP

#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include "cbnetwork/cbnetwork.hpp"

namespace cbnetwork {

class GlobalState {
public:
    std::vector<int> l_values;
    GlobalState(const std::vector<int>& values) : l_values(values) {}
};

class GlobalAttractor {
public:
    std::vector<std::shared_ptr<GlobalState>> l_global_states;
    GlobalAttractor(const std::vector<std::shared_ptr<GlobalState>>& states) : l_global_states(states) {}
};

class GlobalNetwork {
public:
    std::vector<int> l_variables;
    // std::map<int, Function> d_functions; // Future work

    GlobalNetwork() {}

    static void generate_global_network(std::shared_ptr<CBN> o_cbn) {
        if (!o_cbn) return;
        std::cout << "Generating Global Network..." << std::endl;
        std::vector<int> l_variables;
        for (auto& net : o_cbn->l_local_networks) {
            if (!net) continue;
            for (auto& var : net->descriptive_function_variables) {
                if (var) l_variables.push_back(var->index);
            }
        }

        std::cout << "Global variables: [";
        for (size_t i = 0; i < l_variables.size(); ++i) {
            std::cout << l_variables[i] << (i == l_variables.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }

    static void transform_attractor_fields_to_global_states(const std::map<int, std::vector<int>>& l_attractor_fields) {
        // Future work
    }

    static bool test_attractor_fields(std::shared_ptr<CBN> o_cbn) {
        // Future work
        return true;
    }

    static void generate_global_states(const std::vector<int>& o_attractor_field, std::shared_ptr<CBN> o_cbn) {
        // Future work
    }

    bool test_global_dynamic() {
        return true;
    }
};

} // namespace cbnetwork

#endif // GLOBALNETWORK_HPP
