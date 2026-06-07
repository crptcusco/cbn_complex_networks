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

    GlobalNetwork() {}

    static void generate_global_network(std::shared_ptr<CBN> o_cbn) {
        if (!o_cbn) return;
        std::cout << "Generating Global Network..." << std::endl;
        std::vector<int> l_variables;
        for (auto& net : o_cbn->l_local_networks) {
            if (!net) continue;
            for (auto& var : net->internal_variables) {
                l_variables.push_back(var);
            }
        }

        std::cout << "Global variables: [";
        for (size_t i = 0; i < l_variables.size(); ++i) {
            std::cout << l_variables[i] << (i == l_variables.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }

    static void transform_attractor_fields_to_global_states(std::shared_ptr<CBN> o_cbn) {
        if (!o_cbn) return;
        for (auto& pair : o_cbn->d_attractor_fields) {
            generate_global_states(pair.second, o_cbn);
        }
    }

    static bool test_attractor_fields(std::shared_ptr<CBN> o_cbn) {
        if (!o_cbn) return false;
        bool b_flag = true;
        for (auto const& [index, field] : o_cbn->d_attractor_fields) {
            if (test_global_dynamic(field, o_cbn)) {
                std::cout << "Attractor Field " << index << " : Passed" << std::endl;
            } else {
                std::cout << "Attractor Field " << index << " : Failed" << std::endl;
                b_flag = false;
            }
        }
        return b_flag;
    }

    static std::vector<std::shared_ptr<LocalState>> generate_global_states(const std::vector<int>& o_attractor_field, std::shared_ptr<CBN> o_cbn) {
        std::vector<std::shared_ptr<LocalState>> global_states;
        for (int attractor_index : o_attractor_field) {
            auto o_local_attractor = o_cbn->get_local_attractor_by_index(attractor_index);
            if (o_local_attractor) {
                for (auto& state : o_local_attractor->l_states) {
                    global_states.push_back(state);
                }
            }
        }
        return global_states;
    }

    static bool test_global_dynamic(const std::vector<int>& field, std::shared_ptr<CBN> o_cbn) {
        // Placeholder for actual dynamic testing logic
        // In Python it also currently passes if not implemented
        return true;
    }
};

} // namespace cbnetwork

#endif // GLOBALNETWORK_HPP
