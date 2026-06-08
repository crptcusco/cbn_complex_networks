#include "cbnetwork/localnetwork.hpp"
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <omp.h>

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
    external_variables.clear();
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
            else continue;

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

std::vector<int> LocalNetwork::evaluate_all_states_turbo(
    int num_vars,
    const std::vector<std::vector<std::vector<int>>>& cnf_data,
    const std::vector<int>& external_values) {

    long long total_states = 1LL << num_vars;
    std::vector<int> next_states(total_states);

    #pragma omp parallel for
    for (long long state_idx = 0; state_idx < total_states; ++state_idx) {
        int next_state_val = 0;
        for (int i = 0; i < num_vars; ++i) {
            bool var_val = true;
            for (const auto& clause : cnf_data[i]) {
                bool clause_val = false;
                for (int lit : clause) {
                    bool is_neg = (lit < 0);
                    int local_idx = std::abs(lit) - 1;
                    int val;
                    if (local_idx < num_vars) val = (state_idx >> local_idx) & 1;
                    else val = external_values[local_idx - num_vars];

                    if (is_neg) { if (val == 0) { clause_val = true; break; } }
                    else { if (val == 1) { clause_val = true; break; } }
                }
                if (!clause_val) { var_val = false; break; }
            }
            if (var_val) next_state_val |= (1 << i);
        }
        next_states[state_idx] = next_state_val;
    }
    return next_states;
}

void LocalNetwork::find_local_attractors_brute_force_turbo(
    std::shared_ptr<LocalNetwork> local_network,
    const std::vector<std::string>& local_scenes_strings) {

    local_network->local_scenes.clear();
    std::vector<std::string> scenes_to_process = local_scenes_strings;
    if (scenes_to_process.empty()) scenes_to_process.push_back("");

    int num_vars = local_network->internal_variables.size();
    std::map<int, int> var_mapping;
    for (int i = 0; i < num_vars; ++i) var_mapping[local_network->internal_variables[i]] = i;
    for (int i = 0; i < (int)local_network->external_variables.size(); ++i) var_mapping[local_network->external_variables[i]] = num_vars + i;

    std::vector<std::vector<std::vector<int>>> cnf_data(num_vars);
    for (int i = 0; i < num_vars; ++i) {
        for (auto& clause : local_network->descriptive_function_variables[i]->cnf_function) {
            std::vector<int> encoded_clause;
            for (int lit : clause) {
                int encoded_lit = var_mapping[std::abs(lit)] + 1;
                encoded_clause.push_back(lit < 0 ? -encoded_lit : encoded_lit);
            }
            cnf_data[i].push_back(encoded_clause);
        }
    }

    int scene_index = 1;
    for (const auto& scene_str : scenes_to_process) {
        std::vector<int> ext_vals;
        for (size_t i = 0; i < scene_str.length(); ++i) ext_vals.push_back(scene_str[i] - '0');

        auto next_map = evaluate_all_states_turbo(num_vars, cnf_data, ext_vals);

        std::vector<bool> visited(next_map.size(), false);
        std::vector<std::shared_ptr<LocalAttractor>> scene_attractors;

        for (int start_state = 0; start_state < (int)next_map.size(); ++start_state) {
            if (visited[start_state]) continue;

            std::vector<int> path;
            int curr = start_state;
            while (!visited[curr]) {
                visited[curr] = true;
                path.push_back(curr);
                curr = next_map[curr];
            }

            auto it = std::find(path.begin(), path.end(), curr);
            if (it != path.end()) {
                std::vector<std::shared_ptr<LocalState>> l_states;
                for (; it != path.end(); ++it) {
                    std::vector<int> state_vec;
                    for (int bit = 0; bit < num_vars; ++bit) state_vec.push_back((*it >> bit) & 1);
                    l_states.push_back(std::make_shared<LocalState>(state_vec));
                }
                scene_attractors.push_back(std::make_shared<LocalAttractor>(0, scene_attractors.size()+1, l_states, local_network->index, local_network->external_variables, scene_str));
            }
        }

        auto local_scene_obj = std::make_shared<LocalScene>(scene_index++, std::vector<std::string>{scene_str}, local_network->external_variables);
        local_scene_obj->l_attractors = scene_attractors;
        local_network->local_scenes.push_back(local_scene_obj);
    }
    local_network->attractor_count = 0;
    for (auto& s : local_network->local_scenes) local_network->attractor_count += s->l_attractors.size();
}

} // namespace cbnetwork
