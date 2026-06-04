#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/localtemplates.hpp"
#include "cbnetwork/globaltopology.hpp"
#include "cbnetwork/internalvariable.hpp"
#include "cbnetwork/customtext.hpp"
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

void CBN::find_local_attractors_brute_force_turbo() {
    CustomText::make_title("FIND LOCAL ATTRACTORS TURBO");
    #pragma omp parallel for
    for (int i = 0; i < (int)l_local_networks.size(); ++i) {
        auto& net = l_local_networks[i];
        auto scenes = _generate_local_scenes(net);
        LocalNetwork::find_local_attractors_brute_force_turbo(net, scenes);
        process_kind_signal(net);
    }
    _assign_global_indices_to_attractors();
    generate_attractor_dictionary();
    CustomText::make_sub_sub_title("END FIND LOCAL ATTRACTORS TURBO");
}

void CBN::_assign_global_indices_to_attractors() {
    int i_attractor = 1;
    for (auto& net : l_local_networks) {
        for (auto& scene : net->local_scenes) {
            for (auto& attractor : scene->l_attractors) {
                attractor->g_index = i_attractor++;
            }
        }
    }
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
    #pragma omp parallel for
    for (int i = 0; i < (int)l_directed_edges.size(); ++i) {
        auto& edge = l_directed_edges[i];
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
            if (fields[i].empty()) continue;
            for (size_t j = i + 1; j < fields.size(); ++j) {
                if (fields[j].empty()) continue;
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
                    f2.clear(); merged = true; break;
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
    d_attractor_fields.clear();
    for (size_t i = 0; i < current_fields.size(); ++i) d_attractor_fields[i+1] = current_fields[i];
    CustomText::make_sub_sub_title("END MOUNT ATTRACTOR FIELDS TURBO");
}

std::vector<std::vector<int8_t>> CBN::filter_compatible_pairs_kernel(
    const std::vector<std::vector<int>>& fields,
    const std::vector<std::vector<int>>& candidates,
    const std::vector<int>& attr_to_network) {

    size_t n_fields = fields.size();
    size_t n_cands = candidates.size();
    std::vector<std::vector<int8_t>> result(n_fields, std::vector<int8_t>(n_cands, 0));

    #pragma omp parallel for
    for (int i = 0; i < (int)n_fields; ++i) {
        std::set<int> field_nets;
        for (int a : fields[i]) field_nets.insert(attr_to_network[a]);

        for (size_t j = 0; j < n_cands; ++j) {
            int double_check = 0;
            for (int ca : candidates[j]) {
                int cn = attr_to_network[ca];
                if (field_nets.count(cn)) {
                    if (std::find(fields[i].begin(), fields[i].end(), ca) != fields[i].end()) double_check++;
                } else double_check++;
            }
            if (double_check == 2) result[i][j] = 1;
        }
    }
    return result;
}

void CBN::generate_global_scenes() {
    std::vector<int> l_edges_indexes;
    for (const auto& edge : l_directed_edges) l_edges_indexes.push_back(edge->index_variable);
    int n = l_edges_indexes.size();
    l_global_scenes.clear();
    for (int i = 0; i < (1 << n); ++i) {
        std::vector<int> combination;
        for (int bit = 0; bit < n; ++bit) combination.push_back((i >> (n - 1 - bit)) & 1);
        l_global_scenes.push_back(std::make_shared<GlobalScene>(l_edges_indexes, combination));
    }
}

void CBN::count_fields_by_global_scenes() {
    d_global_scenes_count.clear();
    for (auto const& [key, field] : d_attractor_fields) {
        std::map<int, int> d_var_val;
        for (int i_attr : field) {
            auto o_attr = get_local_attractor_by_index(i_attr);
            if (o_attr) {
                for (size_t i = 0; i < o_attr->relation_index.size(); ++i) {
                    d_var_val[o_attr->relation_index[i]] = o_attr->local_scene[i] - '0';
                }
            }
        }
        std::string combo = "";
        for (auto const& [var, val] : d_var_val) combo += std::to_string(val);
        d_global_scenes_count[combo]++;
    }
}

std::shared_ptr<CBN> CBN::cbn_generator(
    int v_topology, int n_local_networks, int n_vars_network,
    int n_input_variables, int n_output_variables,
    int n_max_of_clauses, int n_max_of_literals) {

    auto topo = GlobalTopology::generate_sample_topology(v_topology, n_local_networks);
    LocalNetworkTemplate temp(n_vars_network, n_input_variables, n_output_variables, n_max_of_clauses, n_max_of_literals, v_topology);

    std::vector<std::shared_ptr<LocalNetwork>> networks;
    int var_cont = 1;
    for (int i = 1; i <= n_local_networks; ++i) {
        std::vector<int> ivars;
        for (int j = 0; j < n_vars_network; ++j) ivars.push_back(var_cont++);
        networks.push_back(std::make_shared<LocalNetwork>(i, ivars));
    }

    std::vector<std::shared_ptr<DirectedEdge>> edges;
    int last_var = networks.back()->internal_variables.back() + 1;
    int edge_idx = 1;

    for (auto& rel : topo->l_edges) {
        auto& out_net = networks[rel.first - 1];
        std::vector<int> out_vars;
        for (int p : temp.l_output_var_indexes) out_vars.push_back(out_net->internal_variables[p - 1]);
        edges.push_back(std::make_shared<DirectedEdge>(edge_idx++, last_var++, rel.second, rel.first, out_vars, " OR "));
    }

    for (auto& net : networks) {
        std::vector<std::shared_ptr<DirectedEdge>> inputs;
        for (auto& e : edges) if (e->input_local_network == net->index) inputs.push_back(e);
        net->process_input_signals(inputs);

        for (int ivar : net->internal_variables) {
            int tvar = ivar - ((net->index - 1) * n_vars_network) + n_vars_network;
            auto pre_cnf = temp.d_variable_cnf_function.at(tvar);
            std::vector<std::vector<int>> cnf;
            for (auto& pc : pre_cnf) {
                std::vector<int> clause;
                for (int v : pc) {
                    int lval = std::abs(v) + ((net->index - 3) * n_vars_network) + n_vars_network;
                    if (std::find(net->internal_variables.begin(), net->internal_variables.end(), lval) == net->internal_variables.end()) {
                        if (!net->external_variables.empty()) lval = net->external_variables[0];
                    }
                    clause.push_back(v < 0 ? -lval : lval);
                }
                cnf.push_back(clause);
            }
            net->descriptive_function_variables.push_back(std::make_shared<InternalVariable>(ivar, cnf));
        }
    }
    auto cbn = std::make_shared<CBN>(networks, edges);
    cbn->o_global_topology = topo;
    return cbn;
}

void CBN::show_resume() const {
    std::cout << "CBN Resume: " << l_local_networks.size() << " networks, " << l_directed_edges.size() << " edges." << std::endl;
}

void CBN::show_description() const {
    for (auto& net : l_local_networks) net->show();
}

void CBN::show_global_scenes() const {
    for (auto& s : l_global_scenes) s->show();
}

void CBN::show_stable_attractor_fields() const {
    std::cout << "Attractor Fields: " << d_attractor_fields.size() << std::endl;
}

void CBN::show_local_attractors_dictionary() const {
    for (auto const& [k, v] : d_local_attractors) std::cout << k << " -> (" << std::get<0>(v) << "," << std::get<1>(v) << "," << std::get<2>(v) << ")" << std::endl;
}

void CBN::show_local_attractors() const {
    for (auto& net : l_local_networks) {
        net->show();
        for (auto& scene : net->local_scenes) {
            for (auto& attr : scene->l_attractors) attr->show();
        }
    }
}

void CBN::show_attractor_pairs() const {
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
    for (auto& edge : l_directed_edges) {
        std::cout << "Signal " << edge->index_variable << " Kind: " << edge->kind_signal << std::endl;
    }
}

void CBN::show_stable_attractor_fields_detailed() const {
    show_stable_attractor_fields();
    for (auto const& [id, field] : d_attractor_fields) {
        std::cout << "Field " << id << ": ";
        for (int a : field) std::cout << a << " ";
        std::cout << std::endl;
    }
}

void CBN::show_attractor_fields() const {
    show_stable_attractor_fields();
}

void CBN::show_directed_edges() const {
    for (auto& edge : l_directed_edges) edge->show();
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
