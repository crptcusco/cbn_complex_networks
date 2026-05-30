#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/localtemplates.hpp"
#include "cbnetwork/globaltopology.hpp"
#include "cbnetwork/internalvariable.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <omp.h>
#include "nlohmann/json.hpp"

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

void CBN::find_local_attractors_parallel() {
    #pragma omp parallel for
    for (int i = 0; i < (int)l_local_networks.size(); ++i) {
        auto& net = l_local_networks.at(i);
        if (net) {
            if (net->local_scenes.empty()) {
                auto scenes = _generate_local_scenes(net);
                LocalNetwork::find_local_attractors_brute_force(net, scenes);
            }
            process_kind_signal(net);
        }
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

void CBN::find_compatible_pairs_parallel() {
    #pragma omp parallel for
    for (int i = 0; i < (int)l_directed_edges.size(); ++i) {
        auto& edge = l_directed_edges.at(i);
        if (!edge) continue;
        for (int val : {0, 1}) {
            edge->d_comp_pairs_attractors_by_value[val].clear();
            auto dst_attractors = get_attractors_by_input_signal_value(edge->index_variable, val);

            auto& src_attractors = edge->d_out_value_to_attractor[val];
            edge->d_comp_pairs_attractors_by_value[val].reserve(src_attractors.size() * dst_attractors.size());

            for (auto& src_attr : src_attractors) {
                if (!src_attr) continue;
                for (auto& dst_attr : dst_attractors) {
                    if (!dst_attr) continue;
                    edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_attr->g_index});
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
                                tt_index += std::to_string(state->l_variable_values.at(pos));
                            }
                        }
                    }
                    if (!tt_index.empty() && edge->true_table.count(tt_index)) {
                        signals_in_attractor.insert(edge->true_table.at(tt_index).at(0));
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
                    if (pos < scene_val.length() && scene_val.at(pos) == (signal_value + '0')) {
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

            auto& src_attractors = edge->d_out_value_to_attractor[val];
            edge->d_comp_pairs_attractors_by_value[val].reserve(src_attractors.size() * dst_attractors.size());

            for (auto& src_attr : src_attractors) {
                if (!src_attr) continue;
                for (auto& dst_attr : dst_attractors) {
                    if (!dst_attr) continue;
                    edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_attr->g_index});
                }
            }
        }
    }
}

// OPTIMIZED MEMORY PHASE 3: std::vector<int> instead of std::set<int>
void CBN::mount_stable_attractor_fields() {
    std::vector<std::vector<int>> fields;
    fields.reserve(l_directed_edges.size() * 2);
    for (auto& edge : l_directed_edges) {
        if (!edge) continue;
        for (int val : {0, 1}) {
            for (auto& pair : edge->d_comp_pairs_attractors_by_value[val]) {
                std::vector<int> field = {pair.first, pair.second};
                std::sort(field.begin(), field.end());
                fields.push_back(field);
            }
        }
    }

    if (fields.empty()) return;

    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < fields.size(); ++i) {
            if (fields.at(i).empty()) continue; // Sentinel
            for (size_t j = i + 1; j < fields.size(); ++j) {
                if (fields.at(j).empty()) continue; // Sentinel

                // Efficient sharing check using sorted vectors
                bool share = false;
                auto& f1 = fields.at(i);
                auto& f2 = fields.at(j);

                size_t p1 = 0, p2 = 0;
                while (p1 < f1.size() && p2 < f2.size()) {
                    if (f1[p1] == f2[p2]) { share = true; break; }
                    if (f1[p1] < f2[p2]) p1++;
                    else p2++;
                }

                if (share) {
                    f1.insert(f1.end(), f2.begin(), f2.end());
                    std::sort(f1.begin(), f1.end());
                    f1.erase(std::unique(f1.begin(), f1.end()), f1.end());
                    f2.clear(); // Mark as merged (sentinel)
                    merged = true;
                    break;
                }
            }
            if (merged) break;
        }
    }

    d_attractor_fields.clear();
    int count = 1;
    for (auto& f : fields) {
        if (!f.empty()) {
            d_attractor_fields[count++] = f;
        }
    }
}

void CBN::mount_stable_attractor_fields_parallel() {
    mount_stable_attractor_fields();
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
            std::cout << field.at(i) << (i == field.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
}

std::shared_ptr<CBN> CBN::cbn_generator(
    int v_topology,
    int n_local_networks,
    int n_vars_network,
    int n_input_variables,
    int n_output_variables,
    int n_max_of_clauses,
    int n_max_of_literals) {

    auto o_topo = GlobalTopology::generate_sample_topology(v_topology, n_local_networks);
    LocalNetworkTemplate o_template(n_vars_network, n_input_variables, n_output_variables, n_max_of_clauses, n_max_of_literals, v_topology);

    std::vector<std::shared_ptr<LocalNetwork>> networks;
    int v_cont_var = 1;
    for (int i = 1; i <= n_local_networks; ++i) {
        std::vector<int> internal_vars;
        for (int j = 0; j < n_vars_network; ++j) internal_vars.push_back(v_cont_var++);
        networks.push_back(std::make_shared<LocalNetwork>(i, internal_vars));
    }

    std::vector<std::shared_ptr<DirectedEdge>> edges;
    int i_last_variable = networks.back()->internal_variables.back() + 1;
    int i_directed_edge = 1;

    for (auto& rel : o_topo->l_edges) {
        int out_net_idx = rel.first;
        int in_net_idx = rel.second;

        std::vector<int> out_vars;
        auto& out_net = networks.at(out_net_idx - 1);
        for (int pos : o_template.l_output_var_indexes) {
            out_vars.push_back(out_net->internal_variables.at(pos - 1));
        }

        std::string coup_func = " OR ";
        auto edge = std::make_shared<DirectedEdge>(i_directed_edge++, i_last_variable++, in_net_idx, out_net_idx, out_vars, coup_func);
        edges.push_back(edge);
    }

    for (auto& net : networks) {
        std::vector<std::shared_ptr<DirectedEdge>> inputs;
        for (auto& edge : edges) {
            if (edge->input_local_network == net->index) inputs.push_back(edge);
        }
        net->process_input_signals(inputs);

        for (int i_local_var : net->internal_variables) {
            int i_template_var = i_local_var - ((net->index - 1) * n_vars_network) + n_vars_network;
            auto pre_cnf = o_template.d_variable_cnf_function.at(i_template_var);

            std::vector<std::vector<int>> clauses;
            for (auto& pre_clause : pre_cnf) {
                std::vector<int> clause;
                bool skip = false;
                for (int t_val : pre_clause) {
                    bool pos = t_val > 0;
                    int abs_t = std::abs(t_val);
                    int local_val = abs_t + ((net->index - 3) * n_vars_network) + n_vars_network;

                    if (std::find(net->internal_variables.begin(), net->internal_variables.end(), local_val) == net->internal_variables.end()) {
                        if (net->external_variables.empty()) { skip = true; break; }
                        local_val = net->external_variables.at(0);
                    }
                    clause.push_back(pos ? local_val : -local_val);
                }
                if (!skip && !clause.empty()) clauses.push_back(clause);
            }
            net->descriptive_function_variables.push_back(std::make_shared<InternalVariable>(i_local_var, clauses));
        }
    }

    auto o_cbn = std::make_shared<CBN>(networks, edges);
    o_cbn->o_global_topology = o_topo;
    return o_cbn;
}

void CBN::save_attractor_fields_to_json(const std::string& filepath) {
    using json = nlohmann::json;
    json j_fields = json::array();
    for (auto const& [id, field] : d_attractor_fields) {
        json j_field;
        j_field["field_id"] = id;
        j_field["attractor_indices"] = field;
        j_fields.push_back(j_field);
    }

    std::ofstream out(filepath);
    out << j_fields.dump(4) << std::endl;
}

} // namespace cbnetwork
