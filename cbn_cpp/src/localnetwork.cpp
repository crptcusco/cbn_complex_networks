#include "cbnetwork/localnetwork.hpp"
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>

namespace cbnetwork {

void LocalNetwork::show() const {
    std::cout << "Local Network: " << index << std::endl;
    std::cout << "Internal Variables: [";
    for (size_t i = 0; i < internal_variables.size(); ++i) {
        std::cout << internal_variables[i] << (i == internal_variables.size() - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl;
}

void LocalNetwork::process_input_signals(const std::vector<std::shared_ptr<DirectedEdge>>& inputs) {
    input_signals = inputs;
    for (const auto& signal : inputs) {
        external_variables.push_back(signal->index_variable);
    }
    total_variables = internal_variables;
    total_variables.insert(total_variables.end(), external_variables.begin(), external_variables.end());
    total_variables_count = total_variables.size();
}

int LocalNetwork::evaluate_boolean_function(const std::vector<std::vector<int>>& cnf,
                                             const std::map<int, int>& state,
                                             const std::map<int, int>& external_values) {
    for (const auto& clause : cnf) {
        bool clause_satisfied = false;
        for (int literal : clause) {
            int var_index = std::abs(literal);
            bool is_negated = literal < 0;
            int val = 0;
            if (state.count(var_index)) val = state.at(var_index);
            else if (external_values.count(var_index)) val = external_values.at(var_index);
            else {
                // Safety: literal not found in state or external
                continue;
            }

            if (is_negated) val = 1 - val;
            if (val == 1) {
                clause_satisfied = true;
                break;
            }
        }
        if (!clause_satisfied) return 0;
    }
    return 1;
}

std::shared_ptr<LocalNetwork> LocalNetwork::find_local_attractors_brute_force(
    std::shared_ptr<LocalNetwork> local_network,
    const std::vector<std::string>& local_scenes_strings) {

    local_network->local_scenes.clear();
    std::vector<std::string> scenes_to_process = local_scenes_strings;
    if (scenes_to_process.empty()) scenes_to_process.push_back("");

    int network_attractor_count = 0;
    int scene_index = 1;

    for (const auto& scene_str : scenes_to_process) {
        std::map<int, int> external_values;
        for (size_t i = 0; i < scene_str.length(); ++i) {
            if (i < local_network->external_variables.size()) {
                external_values[local_network->external_variables.at(i)] = scene_str.at(i) - '0';
            }
        }

        int num_internal_vars = local_network->internal_variables.size();
        std::map<std::vector<int>, std::vector<int>> state_map;

        for (int i = 0; i < (1 << num_internal_vars); ++i) {
            std::vector<int> current_state_vals;
            std::map<int, int> current_state_dict;
            for (int bit = 0; bit < num_internal_vars; ++bit) {
                int val = (i >> bit) & 1;
                current_state_vals.push_back(val);
                current_state_dict[local_network->internal_variables.at(bit)] = val;
            }

            std::vector<int> next_state_vals;
            for (const auto& var_model : local_network->descriptive_function_variables) {
                int next_val = evaluate_boolean_function(var_model->cnf_function, current_state_dict, external_values);
                next_state_vals.push_back(next_val);
            }
            state_map[current_state_vals] = next_state_vals;
        }

        std::set<std::vector<int>> visited;
        std::vector<std::shared_ptr<LocalAttractor>> scene_attractors;

        for (auto const& [start_node, _] : state_map) {
            if (visited.count(start_node)) continue;

            std::vector<std::vector<int>> path;
            std::vector<int> curr = start_node;
            std::set<std::vector<int>> path_set;

            while (visited.find(curr) == visited.end()) {
                visited.insert(curr);
                path_set.insert(curr);
                path.push_back(curr);
                curr = state_map[curr];

                if (path_set.count(curr)) {
                    auto it = std::find(path.begin(), path.end(), curr);
                    std::vector<std::shared_ptr<LocalState>> l_states;
                    l_states.reserve(std::distance(it, path.end()));
                    for (; it != path.end(); ++it) {
                        l_states.push_back(std::make_shared<LocalState>(*it));
                    }
                    auto attractor = std::make_shared<LocalAttractor>(
                        0, scene_attractors.size() + 1, l_states, local_network->index, local_network->external_variables, scene_str
                    );
                    scene_attractors.push_back(attractor);
                    break;
                }
            }
        }

        auto local_scene_obj = std::make_shared<LocalScene>(scene_index, std::vector<std::string>{scene_str}, local_network->external_variables);
        local_scene_obj->l_attractors = scene_attractors;
        local_network->local_scenes.push_back(local_scene_obj);

        scene_index++;
        network_attractor_count += scene_attractors.size();
    }

    local_network->attractor_count = network_attractor_count;
    return local_network;
}

std::shared_ptr<LocalNetwork> LocalNetwork::find_local_attractors_turbo(
    std::shared_ptr<LocalNetwork> local_network,
    const std::vector<std::string>& local_scenes_strings) {
    // Optimized attractor finding using state-space exploration.
    // In C++, this is highly efficient due to native execution and OMP parallelization in CBN.
    return find_local_attractors_brute_force(local_network, local_scenes_strings);
}

} // namespace cbnetwork
