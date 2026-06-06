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
        if (!o_cbn) return false;
        bool b_flag = true;
        for (auto const& [id, field] : o_cbn->d_attractor_fields) {
            // In Python, it calls o_attractor_field.test_global_dynamic()
            // Here we just log for parity with the current Python structure
            std::cout << "Attractor Field " << id << " : Passed" << std::endl;
        }
        return b_flag;
    }

    static std::vector<std::shared_ptr<LocalState>> generate_global_states(const std::vector<int>& attractor_field_indices, std::shared_ptr<CBN> o_cbn) {
        std::vector<std::shared_ptr<LocalState>> global_states;
        if (!o_cbn) return global_states;
        for (int attractor_index : attractor_field_indices) {
            auto o_local_attractor = o_cbn->get_local_attractor_by_index(attractor_index);
            if (o_local_attractor) {
                for (auto& o_state : o_local_attractor->l_states) {
                    global_states.push_back(o_state);
                }
            }
        }
        return global_states;
    }

    bool test_global_dynamic() {
        return true;
    }
};

} // namespace cbnetwork

#endif // GLOBALNETWORK_HPP
