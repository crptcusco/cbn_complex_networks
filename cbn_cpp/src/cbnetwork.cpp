#include "cbnetwork/cbnetwork.hpp"
#include <iostream>
#include <algorithm>
#include <set>

namespace cbnetwork {

void CBN::process_output_signals() {
    std::map<int, std::shared_ptr<LocalNetwork>> local_network_dict;
    for (auto& net : l_local_networks) {
        if (net) {
            local_network_dict[net->index] = net;
            net->output_signals.clear();
        }
    }

    for (auto& edge : l_directed_edges) {
        if (!edge) continue;
        int source_network_index = edge->output_local_network;
        if (local_network_dict.count(source_network_index)) {
            local_network_dict[source_network_index]->output_signals.push_back(edge);
        }
    }
}

std::vector<std::string> CBN::_generate_local_scenes(std::shared_ptr<LocalNetwork> o_local_network) {
    if (!o_local_network) return {};
    int external_vars_count = o_local_network->external_variables.size();
    std::vector<std::string> scenes;
    if (external_vars_count > 0) {
        for (int i = 0; i < (1 << external_vars_count); ++i) {
            std::string scene = "";
            for (int bit = 0; bit < external_vars_count; ++bit) {
                scene += ((i >> (external_vars_count - 1 - bit)) & 1) ? "1" : "0";
            }
            scenes.push_back(scene);
        }
        return scenes;
    }
    return {};
}

void CBN::find_local_attractors_sequential() {
    for (auto& net : l_local_networks) {
        if (!net) continue;
        if (net->local_scenes.empty()) {
            auto scenes = _generate_local_scenes(net);
            LocalNetwork::find_local_attractors_brute_force(net, scenes);
        }
    }

    for (auto& net : l_local_networks) {
        if (!net) continue;
        process_kind_signal(net);
    }

    _assign_global_indices_to_attractors();
    generate_attractor_dictionary();
}

void CBN::_assign_global_indices_to_attractors() {
    int i_attractor = 1;
    for (auto& net : l_local_networks) {
        if (!net) continue;
        for (auto& scene : net->local_scenes) {
            if (!scene) continue;
            for (auto& attractor : scene->l_attractors) {
                if (attractor) {
                    attractor->g_index = i_attractor++;
                }
            }
        }
    }
}

void CBN::generate_attractor_dictionary() {
    d_local_attractors.clear();
    for (auto& net : l_local_networks) {
        if (!net) continue;
        for (auto& scene : net->local_scenes) {
            if (!scene) continue;
            for (auto& attractor : scene->l_attractors) {
                if (attractor) {
                    d_local_attractors[attractor->g_index] = {net->index, scene->index, attractor->l_index};
                }
            }
        }
    }
}

void CBN::process_kind_signal(std::shared_ptr<LocalNetwork> o_local_network) {
    if (!o_local_network) return;
    for (auto& edge : l_directed_edges) {
        if (!edge || edge->output_local_network != o_local_network->index) continue;

        edge->d_out_value_to_attractor[0].clear();
        edge->d_out_value_to_attractor[1].clear();

        std::set<int> signals_for_output;

        for (auto& scene : o_local_network->local_scenes) {
            if (!scene) continue;
            std::set<int> signals_in_local_scene;
            for (auto& attractor : scene->l_attractors) {
                if (!attractor) continue;
                std::set<char> signals_in_attractor;
                for (auto& state : attractor->l_states) {
                    if (!state) continue;
                    std::string tt_index = "";
                    for (int out_var : edge->l_output_variables) {
                        auto it = std::find(o_local_network->total_variables.begin(), o_local_network->total_variables.end(), out_var);
                        if (it != o_local_network->total_variables.end()) {
                            int pos = std::distance(o_local_network->total_variables.begin(), it);
                            if (pos < state->l_variable_values.size()) {
                                tt_index += std::to_string(state->l_variable_values[pos]);
                            }
                        }
                    }
                    if (!tt_index.empty() && edge->true_table.count(tt_index)) {
                        signals_in_attractor.insert(edge->true_table.at(tt_index)[0]);
                    }
                }

                if (signals_in_attractor.size() == 1) {
                    int val = *signals_in_attractor.begin() - '0';
                    signals_in_local_scene.insert(val);
                    edge->d_out_value_to_attractor[val].push_back(attractor);
                } else if (!signals_in_attractor.empty()) {
                    signals_in_local_scene.insert(-1);
                }
            }
            for (int s : signals_in_local_scene) signals_for_output.insert(s);
        }

        if (signals_for_output.size() == 1 && *signals_for_output.begin() != -1) {
            edge->kind_signal = 1; // RESTRICTED
        } else if (signals_for_output.count(-1)) {
            edge->kind_signal = 4; // NOT STABLE
        } else if (signals_for_output.size() == 2) {
            edge->kind_signal = 3; // STABLE
        }
    }
}

std::vector<std::shared_ptr<LocalAttractor>> CBN::get_attractors_by_input_signal_value(int index_variable_signal, int signal_value) {
    std::vector<std::shared_ptr<LocalAttractor>> l_attractors;
    for (auto& net : l_local_networks) {
        if (!net) continue;
        for (auto& scene : net->local_scenes) {
            if (!scene) continue;
            auto it = std::find(scene->l_index_signals.begin(), scene->l_index_signals.end(), index_variable_signal);
            if (it != scene->l_index_signals.end()) {
                int pos = std::distance(scene->l_index_signals.begin(), it);
                for (const auto& scene_val : scene->l_values) {
                    if (pos < scene_val.length() && scene_val[pos] == (signal_value + '0')) {
                        for (auto& attr : scene->l_attractors) l_attractors.push_back(attr);
                        break;
                    }
                }
            }
        }
    }
    return l_attractors;
}

void CBN::find_compatible_pairs() {
    for (auto& edge : l_directed_edges) {
        if (!edge) continue;
        for (int val : {0, 1}) {
            edge->d_comp_pairs_attractors_by_value[val].clear();
            auto dst_attractors = get_attractors_by_input_signal_value(edge->index_variable, val);
            for (auto& src_attr : edge->d_out_value_to_attractor[val]) {
                if (!src_attr) continue;
                for (auto& dst_attr : dst_attractors) {
                    if (!dst_attr) continue;
                    edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_attr->g_index});
                }
            }
        }
    }
}

void CBN::mount_stable_attractor_fields() {
    std::vector<std::set<int>> fields;
    for (auto& edge : l_directed_edges) {
        if (!edge) continue;
        for (int val : {0, 1}) {
            for (auto& pair : edge->d_comp_pairs_attractors_by_value[val]) {
                fields.push_back({pair.first, pair.second});
            }
        }
    }

    if (fields.empty()) return;

    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < fields.size(); ++i) {
            for (size_t j = i + 1; j < fields.size(); ++j) {
                bool share = false;
                for (int attr : fields[i]) {
                    if (fields[j].count(attr)) {
                        share = true;
                        break;
                    }
                }
                if (share) {
                    fields[i].insert(fields[j].begin(), fields[j].end());
                    fields.erase(fields.begin() + j);
                    merged = true;
                    break;
                }
            }
            if (merged) break;
        }
    }

    d_attractor_fields.clear();
    for (size_t i = 0; i < fields.size(); ++i) {
        d_attractor_fields[i + 1] = std::vector<int>(fields[i].begin(), fields[i].end());
    }
}

void CBN::show_local_attractors() const {
    for (auto& net : l_local_networks) {
        if (!net) continue;
        for (auto& scene : net->local_scenes) {
            if (!scene) continue;
            for (auto& attr : scene->l_attractors) {
                if (attr) attr->show();
            }
        }
    }
}

void CBN::show_attractor_pairs() const {
    for (auto& edge : l_directed_edges) {
        if (!edge) continue;
        for (int val : {0, 1}) {
            if (edge->d_comp_pairs_attractors_by_value.count(val)) {
                for (auto& pair : edge->d_comp_pairs_attractors_by_value.at(val)) {
                    std::cout << "Pair: (" << pair.first << ", " << pair.second << ")" << std::endl;
                }
            }
        }
    }
}

void CBN::show_stable_attractor_fields() const {
    std::cout << "Number Stable Attractor Fields: " << d_attractor_fields.size() << std::endl;
    for (auto const& [key, field] : d_attractor_fields) {
        std::cout << key << ": [";
        for (size_t i = 0; i < field.size(); ++i) {
            std::cout << field[i] << (i == field.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
}

} // namespace cbnetwork
