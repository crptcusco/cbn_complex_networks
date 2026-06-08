#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/localtemplates.hpp"
#include "cbnetwork/globaltopology.hpp"
#include "cbnetwork/internalvariable.hpp"
#include "cbnetwork/customtext.hpp"
#include "cbnetwork/coupling.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <random>
#include <cmath>
#include <omp.h>
#include <future>
#include "nlohmann/json.hpp"

namespace cbnetwork {

bool CBN::update_network_by_index(std::shared_ptr<LocalNetwork> o_local_network_update) {
    if (!o_local_network_update) return false;
    for (size_t i = 0; i < l_local_networks.size(); ++i) {
        if (l_local_networks[i]->index == o_local_network_update->index) {
            l_local_networks[i] = o_local_network_update;
            return true;
        }
    }
    return false;
}

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

void CBN::order_edges_by_compatibility() {
    if (l_directed_edges.empty()) return;

    auto is_compatible = [](const std::vector<std::shared_ptr<DirectedEdge>>& l_group_base, std::shared_ptr<DirectedEdge> o_group) {
        for (const auto& aux_par : l_group_base) {
            if (aux_par->input_local_network == o_group->input_local_network ||
                aux_par->input_local_network == o_group->output_local_network ||
                aux_par->output_local_network == o_group->output_local_network ||
                aux_par->output_local_network == o_group->input_local_network) {
                return true;
            }
        }
        return false;
    };

    std::vector<std::shared_ptr<DirectedEdge>> l_base = {l_directed_edges[0]};
    std::vector<std::shared_ptr<DirectedEdge>> aux_l_rest_groups(l_directed_edges.begin() + 1, l_directed_edges.end());

    while (!aux_l_rest_groups.empty()) {
        bool found = false;
        for (auto it = aux_l_rest_groups.begin(); it != aux_l_rest_groups.end(); ++it) {
            if (is_compatible(l_base, *it)) {
                l_base.push_back(*it);
                aux_l_rest_groups.erase(it);
                found = true;
                break;
            }
        }
        if (!found) {
            for (auto& e : aux_l_rest_groups) l_base.push_back(e);
            aux_l_rest_groups.clear();
        }
    }
    l_directed_edges = l_base;
}

void CBN::order_edges_by_grade() {
    std::map<int, int> network_degrees;
    for (const auto& net : l_local_networks) network_degrees[net->index] = 0;
    for (const auto& edge : l_directed_edges) {
        network_degrees[edge->input_local_network]++;
        network_degrees[edge->output_local_network]++;
    }

    auto calculate_edge_grade = [&](std::shared_ptr<DirectedEdge> edge) {
        return network_degrees[edge->input_local_network] + network_degrees[edge->output_local_network];
    };

    std::sort(l_directed_edges.begin(), l_directed_edges.end(), [&](std::shared_ptr<DirectedEdge> a, std::shared_ptr<DirectedEdge> b) {
        return calculate_edge_grade(a) > calculate_edge_grade(b);
    });

    auto is_adjacent = [](std::shared_ptr<DirectedEdge> edge1, std::shared_ptr<DirectedEdge> edge2) {
        return edge1->input_local_network == edge2->input_local_network ||
               edge1->input_local_network == edge2->output_local_network ||
               edge1->output_local_network == edge2->output_local_network ||
               edge1->output_local_network == edge2->input_local_network;
    };

    if (l_directed_edges.empty()) return;

    std::vector<std::shared_ptr<DirectedEdge>> ordered_edges = {l_directed_edges[0]};
    l_directed_edges.erase(l_directed_edges.begin());

    while (!l_directed_edges.empty()) {
        bool found = false;
        for (auto it = l_directed_edges.begin(); it != l_directed_edges.end(); ++it) {
            if (is_adjacent(ordered_edges.back(), *it)) {
                ordered_edges.push_back(*it);
                l_directed_edges.erase(it);
                found = true;
                break;
            }
        }
        if (!found) {
            ordered_edges.push_back(l_directed_edges[0]);
            l_directed_edges.erase(l_directed_edges.begin());
        }
    }
    l_directed_edges = ordered_edges;
}

void CBN::disorder_edges() {
    if (l_directed_edges.size() < 2) return;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(l_directed_edges.begin(), l_directed_edges.end(), g);

    auto have_common_vertex = [](std::shared_ptr<DirectedEdge> edge1, std::shared_ptr<DirectedEdge> edge2) {
        return edge1->input_local_network == edge2->input_local_network ||
               edge1->input_local_network == edge2->output_local_network ||
               edge1->output_local_network == edge2->output_local_network ||
               edge1->output_local_network == edge2->input_local_network;
    };

    if (have_common_vertex(l_directed_edges[0], l_directed_edges[1])) {
        for (size_t i = 2; i < l_directed_edges.size(); ++i) {
            if (!have_common_vertex(l_directed_edges[0], l_directed_edges[i])) {
                std::swap(l_directed_edges[1], l_directed_edges[i]);
                break;
            }
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
    } else {
        scenes.push_back("");
    }
    return scenes;
}

void CBN::find_local_attractors_sequential() {
    CustomText::make_title("FIND LOCAL ATTRACTORS");
    for (auto& net : l_local_networks) {
        auto scenes = _generate_local_scenes(net);
        LocalNetwork::find_local_attractors_brute_force(net, scenes);
        process_kind_signal(net);
    }
    _assign_global_indices_to_attractors();
    generate_attractor_dictionary();
    CustomText::make_sub_sub_title("END FIND LOCAL ATTRACTORS");
}

void CBN::find_local_attractors_parallel() {
    CustomText::make_title("FIND LOCAL ATTRACTORS PARALLEL");
    #pragma omp parallel for
    for (int i = 0; i < (int)l_local_networks.size(); ++i) {
        auto& net = l_local_networks[i];
        auto scenes = _generate_local_scenes(net);
        LocalNetwork::find_local_attractors_brute_force(net, scenes);
        process_kind_signal(net);
    }
    _assign_global_indices_to_attractors();
    generate_attractor_dictionary();
    CustomText::make_sub_sub_title("END FIND LOCAL ATTRACTORS PARALLEL");
}

void CBN::find_local_attractors_parallel_with_weights() {
    CustomText::make_title("FIND LOCAL ATTRACTORS WEIGHTED");
    struct WeightedNet { long weight; std::shared_ptr<LocalNetwork> net; };
    std::vector<WeightedNet> weighted_nets;
    for (auto& net : l_local_networks) {
        long w = (long)net->internal_variables.size() * (1L << net->external_variables.size());
        weighted_nets.push_back({w, net});
    }
    std::sort(weighted_nets.begin(), weighted_nets.end(), [](const auto& a, const auto& b){ return a.weight > b.weight; });

    int num_threads = omp_get_max_threads();
    std::vector<std::vector<std::shared_ptr<LocalNetwork>>> buckets(num_threads);
    std::vector<long> bucket_weights(num_threads, 0);
    for (auto& wn : weighted_nets) {
        int best_idx = 0; long min_w = bucket_weights[0];
        for(int i=1; i<num_threads; ++i) { if(bucket_weights[i] < min_w) { min_w=bucket_weights[i]; best_idx=i; } }
        buckets[best_idx].push_back(wn.net);
        bucket_weights[best_idx] += wn.weight;
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (auto& net : buckets[tid]) {
            auto scenes = _generate_local_scenes(net);
            LocalNetwork::find_local_attractors_brute_force(net, scenes);
            process_kind_signal(net);
        }
    }
    _assign_global_indices_to_attractors();
    generate_attractor_dictionary();
    CustomText::make_sub_sub_title("END FIND LOCAL ATTRACTORS WEIGHTED");
}

void CBN::generate_attractor_dictionary() {
    d_local_attractors.clear();
    d_local_attractors_ptr.clear();
    for (auto& net : l_local_networks) {
        for (auto& scene : net->local_scenes) {
            for (auto& attractor : scene->l_attractors) {
                d_local_attractors[attractor->g_index] = {net->index, scene->index, attractor->l_index};
                d_local_attractors_ptr[attractor->g_index] = attractor;
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
            std::set<int> signals_in_local_scene;
            for (auto& attractor : scene->l_attractors) {
                std::set<char> signals_in_attractor;
                for (auto& state : attractor->l_states) {
                    std::string tt_index = "";
                    for (int out_var : edge->l_output_variables) {
                        auto it = std::find(o_local_network->total_variables.begin(), o_local_network->total_variables.end(), out_var);
                        if (it != o_local_network->total_variables.end()) {
                            int pos = std::distance(o_local_network->total_variables.begin(), it);
                            tt_index += std::to_string(state->l_variable_values.at(pos));
                        }
                    }
                    if (edge->true_table.count(tt_index)) signals_in_attractor.insert(edge->true_table.at(tt_index).at(0));
                }

                if (signals_in_attractor.size() == 1) {
                    int val = *signals_in_attractor.begin() - '0';
                    signals_in_local_scene.insert(val);
                    edge->d_out_value_to_attractor[val].push_back(attractor);
                } else {
                    signals_in_local_scene.insert(-1);
                }
            }
            for (int s : signals_in_local_scene) signals_for_output.insert(s);
        }

        if (signals_for_output.size() == 1 && *signals_for_output.begin() != -1) edge->kind_signal = 1;
        else if (signals_for_output.count(-1)) edge->kind_signal = 4;
        else if (signals_for_output.size() == 2) edge->kind_signal = 3;
    }
}

std::vector<std::shared_ptr<LocalAttractor>> CBN::get_attractors_by_input_signal_value(int index_variable_signal, int signal_value) {
    std::vector<std::shared_ptr<LocalAttractor>> l_attractors;
    for (auto& net : l_local_networks) {
        for (auto& scene : net->local_scenes) {
            auto it = std::find(scene->l_index_signals.begin(), scene->l_index_signals.end(), index_variable_signal);
            if (it != scene->l_index_signals.end()) {
                int pos = std::distance(scene->l_index_signals.begin(), it);
                if (scene->l_values.at(0).at(pos) == (signal_value + '0')) {
                    for (auto& attr : scene->l_attractors) l_attractors.push_back(attr);
                }
            }
        }
    }
    return l_attractors;
}

std::shared_ptr<LocalAttractor> CBN::get_local_attractor_by_index(int i_attractor) {
    if (d_local_attractors_ptr.count(i_attractor)) return d_local_attractors_ptr.at(i_attractor);
    return nullptr;
}

std::shared_ptr<LocalNetwork> CBN::get_network_by_index(int index) {
    for (auto& net : l_local_networks) if (net->index == index) return net;
    return nullptr;
}

void CBN::find_compatible_pairs() {
    for (auto& edge : l_directed_edges) {
        for (int val : {0, 1}) {
            edge->d_comp_pairs_attractors_by_value[val].clear();
            auto dst_attractors = get_attractors_by_input_signal_value(edge->index_variable, val);
            for (auto& src_attr : edge->d_out_value_to_attractor[val]) {
                for (auto& dst_attr : dst_attractors) {
                    edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_attr->g_index});
                }
            }
        }
    }
}

void CBN::find_compatible_pairs_parallel() {
    CustomText::make_title("FIND COMPATIBLE ATTRACTOR PAIRS (PARALLEL)");
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
    CustomText::make_sub_sub_title("END FIND COMPATIBLE ATTRACTOR PAIRS (PARALLEL)");
}

void CBN::find_compatible_pairs_parallel_with_weights() {
    CustomText::make_title("FIND COMPATIBLE PAIRS WEIGHTED");
    struct WeightedEdge { long weight; std::shared_ptr<DirectedEdge> edge; };
    std::vector<WeightedEdge> weighted_edges;
    for (auto& edge : l_directed_edges) {
        long w = (long)(edge->d_out_value_to_attractor[0].size() + edge->d_out_value_to_attractor[1].size());
        weighted_edges.push_back({w, edge});
    }
    std::sort(weighted_edges.begin(), weighted_edges.end(), [](const auto& a, const auto& b){ return a.weight > b.weight; });

    int num_threads = omp_get_max_threads();
    std::vector<std::vector<std::shared_ptr<DirectedEdge>>> buckets(num_threads);
    std::vector<long> bucket_weights(num_threads, 0);
    for (auto& we : weighted_edges) {
        int best_idx = 0; long min_w = bucket_weights[0];
        for(int i=1; i<num_threads; ++i) { if(bucket_weights[i] < min_w) { min_w=bucket_weights[i]; best_idx=i; } }
        buckets[best_idx].push_back(we.edge);
        bucket_weights[best_idx] += we.weight;
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (auto& edge : buckets[tid]) {
            for (int val : {0, 1}) {
                edge->d_comp_pairs_attractors_by_value[val].clear();
                auto dst_attractors = get_attractors_by_input_signal_value(edge->index_variable, val);
                for (auto& src_attr : edge->d_out_value_to_attractor[val]) {
                    for (auto& dst_attr : dst_attractors) {
                        edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_attr->g_index});
                    }
                }
            }
        }
    }
    CustomText::make_sub_sub_title("END FIND COMPATIBLE PAIRS WEIGHTED");
}

void CBN::find_compatible_pairs_turbo() {
    CustomText::make_title("FIND COMPATIBLE PAIRS TURBO");
    find_compatible_pairs_parallel();
    CustomText::make_sub_sub_title("END FIND COMPATIBLE PAIRS TURBO");
}

void CBN::mount_stable_attractor_fields_parallel() {
    mount_stable_attractor_fields();
}

void CBN::mount_stable_attractor_fields_parallel_chunks() {
    CustomText::make_title("MOUNT ATTRACTOR FIELDS CHUNKS");
    mount_stable_attractor_fields();
    CustomText::make_sub_sub_title("END MOUNT ATTRACTOR FIELDS CHUNKS");
}

void CBN::mount_stable_attractor_fields() {
    std::vector<std::vector<int>> fields;
    for (auto& edge : l_directed_edges) {
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
            if (fields.at(i).empty()) continue;
            for (size_t j = i + 1; j < fields.size(); ++j) {
                if (fields.at(j).empty()) continue;

                bool share = false;
                auto& f1 = fields[i]; auto& f2 = fields[j];
                size_t p1 = 0, p2 = 0;
                while (p1 < f1.size() && p2 < f2.size()) {
                    if (f1[p1] == f2[p2]) { share = true; break; }
                    if (f1[p1] < f2[p2]) p1++; else p2++;
                }
                if (share) {
                    f1.insert(f1.end(), f2.begin(), f2.end());
                    std::sort(f1.begin(), f1.end());
                    f1.erase(std::unique(f1.begin(), f1.end()), f1.end());
                    f2.clear();
                    merged = true;
                    break;
                }
            }
            if (merged) break;
        }
    }
    d_attractor_fields.clear();
    int count = 1;
    for (auto& f : fields) if (!f.empty()) d_attractor_fields[count++] = f;
}

void CBN::mount_stable_attractor_fields_turbo() {
    CustomText::make_title("MOUNT ATTRACTOR FIELDS TURBO");
    order_edges_by_compatibility();

    int max_attr_idx = 0;
    for (auto const& [idx, _] : d_local_attractors) max_attr_idx = std::max(max_attr_idx, idx);
    std::vector<int> attr_to_network(max_attr_idx + 1);
    for (auto const& [idx, info] : d_local_attractors) attr_to_network[idx] = std::get<0>(info);

    std::vector<std::vector<int>> current_fields;
    if (l_directed_edges.empty()) return;

    auto first_edge = l_directed_edges[0];
    for (int val : {0, 1}) {
        for (auto& p : first_edge->d_comp_pairs_attractors_by_value[val]) current_fields.push_back({p.first, p.second});
    }

    for (size_t i = 1; i < l_directed_edges.size(); ++i) {
        std::vector<std::vector<int>> candidates;
        for (int val : {0, 1}) {
            for (auto& p : l_directed_edges[i]->d_comp_pairs_attractors_by_value[val]) candidates.push_back({p.first, p.second});
        }
        if (candidates.empty() || current_fields.empty()) { current_fields.clear(); break; }

        auto comp_matrix = filter_compatible_pairs_kernel(current_fields, candidates, attr_to_network);
        std::vector<std::vector<int>> next_fields;
        for (size_t f_idx = 0; f_idx < current_fields.size(); ++f_idx) {
            for (size_t c_idx = 0; c_idx < candidates.size(); ++c_idx) {
                if (comp_matrix[f_idx][c_idx]) {
                    std::vector<int> new_f = current_fields[f_idx];
                    new_f.insert(new_f.end(), candidates[c_idx].begin(), candidates[c_idx].end());
                    std::sort(new_f.begin(), new_f.end());
                    new_f.erase(std::unique(new_f.begin(), new_f.end()), new_f.end());
                    next_fields.push_back(new_f);
                }
            }
        }
        current_fields = next_fields;
    }
    int total_f_pairs = 0;
    for (const auto& f : current_fields) total_f_pairs += f.size();
    std::cout << "Total compatible attractor pairs in fields: " << total_f_pairs << std::endl;
}

std::vector<std::vector<int8_t>> CBN::filter_compatible_pairs_kernel(
    const std::vector<std::vector<int>>& fields,
    const std::vector<std::vector<int>>& candidates,
    const std::vector<int>& attr_to_network) {

    size_t n_fields = fields.size();
    size_t n_cands = candidates.size();
    std::vector<std::vector<int8_t>> result(n_fields, std::vector<int8_t>(n_cands, 0));

    int max_net = 0;
    for (int net_idx : attr_to_network) if (net_idx > max_net) max_net = net_idx;

    #pragma omp parallel for
    for (int i = 0; i < (int)n_fields; ++i) {
        std::vector<int> net_to_attr(max_net + 1, -1);
        for (int attr : fields[i]) {
            net_to_attr[attr_to_network[attr]] = attr;
        }
        for (int j = 0; j < (int)n_cands; ++j) {
            bool compatible = true;
            for (int attr : candidates[j]) {
                int net = attr_to_network[attr];
                if (net_to_attr[net] != -1 && net_to_attr[net] != attr) {
                    compatible = false;
                    break;
                }
            }
            if (compatible) result[i][j] = 1;
        }
    }
    return result;
}

std::shared_ptr<CBN> CBN::cbn_generator(
    int v_topology,
    int n_local_networks,
    int n_vars_network,
    int n_input_variables,
    int n_output_variables,
    int n_max_of_clauses,
    int n_max_of_literals,
    int n_edges,
    std::shared_ptr<CouplingStrategy> coupling_strategy) {

    auto o_topo = GlobalTopology::generate_sample_topology(v_topology, n_local_networks, n_edges);
    LocalNetworkTemplate o_template(n_vars_network, n_input_variables, n_output_variables, n_max_of_clauses, n_max_of_literals, v_topology);

    std::vector<std::shared_ptr<LocalNetwork>> networks;
    int var_cont = 1;
    for (int i = 1; i <= n_local_networks; ++i) {
        std::vector<int> ivars;
        for (int j = 0; j < n_vars_network; ++j) ivars.push_back(var_cont++);
        networks.push_back(std::make_shared<LocalNetwork>(i, ivars));
    }

    std::vector<std::shared_ptr<DirectedEdge>> edges;
    int i_last_variable = var_cont;
    int i_directed_edge = 1;

    for (auto& rel : o_topo->l_edges) {
        int out_net_idx = rel.first;
        int in_net_idx = rel.second;

        auto& out_net = networks.at(out_net_idx - 1);
        std::vector<int> out_vars;
        for (int pos : o_template.l_output_var_indexes) {
            out_vars.push_back(out_net->internal_variables.at(pos - 1));
        }

        std::string coup_func = coupling_strategy->generate_coupling_function(out_vars);
        auto edge = std::make_shared<DirectedEdge>(i_directed_edge++, i_last_variable, in_net_idx, out_net_idx, out_vars, coup_func);

        auto coupling_cnf = coupling_strategy->to_cnf(out_vars, i_last_variable);
        auto coupling_variable = std::make_shared<InternalVariable>(i_last_variable, coupling_cnf);

        if (std::find(out_net->internal_variables.begin(), out_net->internal_variables.end(), i_last_variable) == out_net->internal_variables.end()) {
            out_net->internal_variables.push_back(i_last_variable);
        }
        out_net->descriptive_function_variables.push_back(coupling_variable);

        i_last_variable++;
        edges.push_back(edge);
    }

    for (auto& net : networks) {
        auto inputs = find_input_edges_by_network_index(net->index, edges);
        net->process_input_signals(inputs);

        for (int ivar : net->internal_variables) {
            int tvar = ivar - ((net->index - 1) * n_vars_network);
            if (tvar < 1 || tvar > n_vars_network) continue;

            auto it_cnf = o_template.d_variable_cnf_function.find(tvar + n_vars_network);
            if (it_cnf == o_template.d_variable_cnf_function.end()) continue;

            const auto& pre_cnf = it_cnf->second;
            std::vector<std::vector<int>> cnf;
            for (auto& pc : pre_cnf) {
                std::vector<int> clause;
                for (int v : pc) {
                    int abs_v = std::abs(v);
                    int lval;
                    if (abs_v >= n_vars_network + 1 && abs_v <= n_vars_network * 2) {
                        lval = (abs_v - n_vars_network) + (net->index - 1) * n_vars_network;
                    } else if (abs_v > n_vars_network * 2) {
                        int signal_idx = abs_v - (n_vars_network * 2) - 1;
                        if (signal_idx < (int)net->input_signals.size()) {
                            lval = net->input_signals[signal_idx]->index_variable;
                        } else {
                            lval = net->internal_variables[0];
                        }
                    } else {
                        lval = abs_v;
                    }
                    clause.push_back(v < 0 ? -lval : lval);
                }
                cnf.push_back(clause);
            }
            net->descriptive_function_variables.push_back(std::make_shared<InternalVariable>(ivar, cnf));
        }
    }
    auto cbn = std::make_shared<CBN>(networks, edges);
    cbn->o_global_topology = o_topo;
    return cbn;
}

void CBN::show_resume() const {
    CustomText::make_title("CBN Detailed Resume");
    CustomText::make_sub_sub_title("Main Characteristics");
    CustomText::print_simple_line();
    std::cout << "Number of local networks: " << l_local_networks.size() << std::endl;
    std::cout << "Topology Type: " << (o_global_topology ? std::to_string(o_global_topology->v_topology) : "N/A") << std::endl;

    CustomText::make_sub_sub_title("Indicators");
    CustomText::print_simple_line();
    std::cout << "Number of local attractors: " << get_n_local_attractors() << std::endl;
    std::cout << "Number of attractor fields: " << get_n_attractor_fields() << std::endl;
    CustomText::print_simple_line();
}

void CBN::show_description() const {
    CustomText::make_title("CBN description");
    std::string nets = "Local Networks: [";
    for(size_t i=0; i<l_local_networks.size(); ++i) nets += std::to_string(l_local_networks[i]->index) + (i==l_local_networks.size()-1?"":", ");
    nets += "]";
    CustomText::make_sub_title(nets);
    for (auto& net : l_local_networks) if(net) net->show();
}

void CBN::show_global_scenes() const {
    CustomText::make_sub_title("LIST OF GLOBAL SCENES");
    for (auto& scene : l_global_scenes) if(scene) scene->show();
}

void CBN::show_stable_attractor_fields() const {
    CustomText::print_duplex_line();
    std::cout << "Show the list of attractor fields" << std::endl;
    std::cout << "Number Stable Attractor Fields: " << d_attractor_fields.size() << std::endl;
    for (auto const& [key, field] : d_attractor_fields) {
        CustomText::print_simple_line();
        std::cout << key << ": [";
        for (size_t i = 0; i < field.size(); ++i) {
            std::cout << field.at(i) << (i == field.size() - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
}

void CBN::show_local_attractors_dictionary() const {
    CustomText::make_title("Global Dictionary of Local Attractors");
    for (auto const& [key, val] : d_local_attractors) {
        std::cout << key << " -> (" << std::get<0>(val) << ", " << std::get<1>(val) << ", " << std::get<2>(val) << ")" << std::endl;
    }
}

void CBN::show_local_attractors() const {
    CustomText::make_title("SHOW LOCAL ATTRACTORS");
    for (auto& net : l_local_networks) {
        net->show();
        for (auto& scene : net->local_scenes) {
            for (auto& attr : scene->l_attractors) attr->show();
        }
    }
}

void CBN::show_attractor_pairs() const {
    CustomText::make_title("SHOW ATTRACTOR PAIRS");
    for (auto& edge : l_directed_edges) {
        std::cout << "Edge " << edge->index << ": " << edge->output_local_network << " -> " << edge->input_local_network << std::endl;
        for (int val : {0, 1}) {
            if (edge->d_comp_pairs_attractors_by_value.count(val)) {
                for (auto& p : edge->d_comp_pairs_attractors_by_value.at(val)) {
                    std::cout << "  (" << p.first << ", " << p.second << ")" << std::endl;
                }
            }
        }
    }
}

void CBN::show_coupled_signals_kind() const {
    CustomText::print_duplex_line();
    std::cout << "SHOW THE COUPLED SIGNALS KINDS" << std::endl;
    int n_restricted = 0;
    for (auto& edge : l_directed_edges) {
        if (!edge) continue;
        std::cout << "SIGNAL: " << edge->index_variable
                  << ", RELATION: " << edge->output_local_network << " -> " << edge->input_local_network
                  << ", KIND: " << edge->kind_signal << std::endl;
        if (edge->kind_signal == 1) n_restricted++;
    }
    std::cout << "Number of restricted signals: " << n_restricted << std::endl;
}

void CBN::show_stable_attractor_fields_detailed() const {
    show_stable_attractor_fields();
    for (auto const& [id, field] : d_attractor_fields) {
        std::cout << "Field " << id << ": ";
        for (int a : field) {
            if (d_local_attractors.count(a)) {
                auto val = d_local_attractors.at(a);
                std::cout << a << " -> (" << std::get<0>(val) << ", " << std::get<1>(val) << ", " << std::get<2>(val) << ") ";
            }
        }
        std::cout << std::endl;
    }
}

void CBN::show_attractor_fields() const {
    show_stable_attractor_fields();
}

void CBN::show_directed_edges() const {
    CustomText::print_duplex_line();
    std::cout << "SHOW THE DIRECTED EDGES OF THE CBN" << std::endl;
    for (auto& edge : l_directed_edges) {
        if (edge) edge->show();
    }
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

void CBN::save_network_to_json(const std::string& filepath) const {
    using json = nlohmann::json;
    json j_cbn;

    json j_nets = json::array();
    for (const auto& net : l_local_networks) {
        json j_net;
        j_net["index"] = net->index;
        j_net["internal_variables"] = net->internal_variables;

        json j_vars = json::array();
        for (const auto& var : net->descriptive_function_variables) {
            json j_var;
            j_var["index"] = var->index;
            j_var["cnf"] = var->cnf_function;
            j_vars.push_back(j_var);
        }
        j_net["logic"] = j_vars;
        j_nets.push_back(j_net);
    }
    j_cbn["local_networks"] = j_nets;

    json j_edges = json::array();
    for (const auto& edge : l_directed_edges) {
        json j_edge;
        j_edge["index"] = edge->index;
        j_edge["signal_var"] = edge->index_variable;
        j_edge["source"] = edge->output_local_network;
        j_edge["target"] = edge->input_local_network;
        j_edge["output_vars"] = edge->l_output_variables;
        j_edge["function"] = edge->coupling_function;
        j_edges.push_back(j_edge);
    }
    j_cbn["directed_edges"] = j_edges;

    std::ofstream out(filepath);
    out << j_cbn.dump(4) << std::endl;
}

double CBN::calculate_canalization_density() const {
    int total_internal_vars = 0;
    std::set<int> used_vars;

    for (const auto& net : l_local_networks) {
        total_internal_vars += net->internal_variables.size();
        for (const auto& var : net->descriptive_function_variables) {
            for (const auto& clause : var->cnf_function) {
                for (int literal : clause) {
                    used_vars.insert(std::abs(literal));
                }
            }
        }
    }

    if (total_internal_vars == 0) return 0.0;

    int redundant_vars = 0;
    for (const auto& net : l_local_networks) {
        for (int ivar : net->internal_variables) {
            if (used_vars.find(ivar) == used_vars.end()) {
                redundant_vars++;
            }
        }
    }

    return (double)redundant_vars / total_internal_vars;
}

std::vector<std::shared_ptr<DirectedEdge>> CBN::find_output_edges_by_network_index(int index, const std::vector<std::shared_ptr<DirectedEdge>>& edges) {
    std::vector<std::shared_ptr<DirectedEdge>> result;
    for(auto& e : edges) if(e->output_local_network == index) result.push_back(e);
    return result;
}

std::vector<std::shared_ptr<DirectedEdge>> CBN::find_input_edges_by_network_index(int index, const std::vector<std::shared_ptr<DirectedEdge>>& edges) {
    std::vector<std::shared_ptr<DirectedEdge>> result;
    for(auto& e : edges) if(e->input_local_network == index) result.push_back(e);
    return result;
}

int CBN::get_n_local_attractors() const {
    int count = 0;
    for (const auto& net : l_local_networks) count += net->attractor_count;
    return count;
}

int CBN::get_n_attractor_fields() const {
    return d_attractor_fields.size();
}

void CBN::_assign_global_indices_to_attractors() {
    int g_index = 1;
    for (auto& net : l_local_networks) {
        for (auto& scene : net->local_scenes) {
            for (auto& attr : scene->l_attractors) {
                attr->g_index = g_index++;
            }
        }
    }
}

void CBN::generate_global_scenes() {}
void CBN::count_fields_by_global_scenes() {}

} // namespace cbnetwork
