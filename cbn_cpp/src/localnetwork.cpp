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

    local_network->local_scenes.clear();
    std::vector<std::string> scenes_to_process = local_scenes_strings;
    if (scenes_to_process.empty()) scenes_to_process.push_back("");

    int num_vars = local_network->internal_variables.size();
    int total_states = 1 << num_vars;

    // Pre-calculate CNF structure for turbo speed
    struct Literal { int local_idx; bool is_neg; bool is_external; };
    struct Clause { std::vector<Literal> literals; };
    struct VarFunc { std::vector<Clause> clauses; };

    std::vector<VarFunc> var_functions;
    std::map<int, int> global_to_local;
    for(int i=0; i<num_vars; ++i) global_to_local[local_network->internal_variables[i]] = i;

    std::map<int, int> global_to_ext;
    for(int i=0; i<(int)local_network->external_variables.size(); ++i) global_to_ext[local_network->external_variables[i]] = i;

    for (const auto& var_model : local_network->descriptive_function_variables) {
        VarFunc vf;
        for (const auto& clause_cnf : var_model->cnf_function) {
            Clause c;
            for (int lit : clause_cnf) {
                Literal l;
                int g_idx = std::abs(lit);
                l.is_neg = lit < 0;
                if (global_to_local.count(g_idx)) {
                    l.local_idx = global_to_local[g_idx];
                    l.is_external = false;
                } else {
                    l.local_idx = global_to_ext[g_idx];
                    l.is_external = true;
                }
                c.literals.push_back(l);
            }
            vf.clauses.push_back(c);
        }
        var_functions.push_back(vf);
    }

    int scene_index = 1;
    int network_attractor_count = 0;

    for (const auto& scene_str : scenes_to_process) {
        std::vector<int> external_vals(local_network->external_variables.size(), 0);
        for (size_t i = 0; i < scene_str.length() && i < external_vals.size(); ++i) {
            external_vals[i] = scene_str[i] - '0';
        }

        std::vector<int> next_states(total_states);

        #pragma omp parallel for
        for (int s = 0; s < total_states; ++s) {
            int next_s = 0;
            for (int v = 0; v < num_vars; ++v) {
                bool var_val = true;
                for (const auto& clause : var_functions[v].clauses) {
                    bool clause_val = false;
                    for (const auto& lit : clause.literals) {
                        int val = lit.is_external ? external_vals[lit.local_idx] : (s >> lit.local_idx) & 1;
                        if (lit.is_neg) val = 1 - val;
                        if (val == 1) { clause_val = true; break; }
                    }
                    if (!clause_val) { var_val = false; break; }
                }
                if (var_val) next_s |= (1 << v);
            }
            next_states[s] = next_s;
        }

        std::vector<bool> visited(total_states, false);
        std::vector<std::shared_ptr<LocalAttractor>> scene_attractors;

        for (int i = 0; i < total_states; ++i) {
            if (visited[i]) continue;
            std::vector<int> path;
            int curr = i;
            while (!visited[curr]) {
                visited[curr] = true;
                path.push_back(curr);
                curr = next_states[curr];
            }

            auto it = std::find(path.begin(), path.end(), curr);
            if (it != path.end()) {
                std::vector<std::shared_ptr<LocalState>> l_states;
                for (; it != path.end(); ++it) {
                    std::vector<int> state_vals(num_vars);
                    for(int b=0; b<num_vars; ++b) state_vals[b] = (*it >> b) & 1;
                    l_states.push_back(std::make_shared<LocalState>(state_vals));
                }
                auto attractor = std::make_shared<LocalAttractor>(
                    0, scene_attractors.size() + 1, l_states, local_network->index, local_network->external_variables, scene_str
                );
                scene_attractors.push_back(attractor);
            }
        }

        auto local_scene_obj = std::make_shared<LocalScene>(scene_index++, std::vector<std::string>{scene_str}, local_network->external_variables);
        local_scene_obj->l_attractors = scene_attractors;
        local_network->local_scenes.push_back(local_scene_obj);
        network_attractor_count += scene_attractors.size();
    }
    local_network->attractor_count = network_attractor_count;
    return local_network;
}
}
