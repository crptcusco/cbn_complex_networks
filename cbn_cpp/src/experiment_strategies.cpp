#include "cbnetwork/experiment_strategies.hpp"
#include <omp.h>
#include <chrono>
#include <sys/resource.h>
#include <future>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <functional>

namespace cbnetwork {

using namespace std::chrono;

static long get_max_rss() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

// Memory-efficient vector merge
static std::vector<int> merge_two_sorted_vectors(const std::vector<int>& v1, const std::vector<int>& v2) {
    std::vector<int> res;
    res.reserve(v1.size() + v2.size());
    std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(res));
    return res;
}

static std::vector<std::vector<int>> merge_field_lists_sequentially(std::vector<std::vector<int>> fields) {
    if (fields.empty()) return {};
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < fields.size(); ++i) {
            if (fields.at(i).empty()) continue;
            for (size_t j = i + 1; j < fields.size(); ++j) {
                if (fields.at(j).empty()) continue;

                auto& f1 = fields.at(i);
                auto& f2 = fields.at(j);
                bool share = false;
                size_t p1 = 0, p2 = 0;
                while (p1 < f1.size() && p2 < f2.size()) {
                    if (f1[p1] == f2[p2]) { share = true; break; }
                    if (f1[p1] < f2[p2]) p1++; else p2++;
                }

                if (share) {
                    f1.insert(f1.end(), f2.begin(), f2.end());
                    std::sort(f1.begin(), f1.end());
                    f1.erase(std::unique(f1.begin(), f1.end()), f1.end());
                    f2.clear(); // Sentinel
                    merged = true;
                    break;
                }
            }
            if (merged) break;
        }
    }
    std::vector<std::vector<int>> final_res;
    for (auto& f : fields) if (!f.empty()) final_res.push_back(f);
    return final_res;
}

// 1. TRADITIONAL EXPERIMENT
ExperimentResults TraditionalExperiment::run(std::shared_ptr<CBN> cbn) {
    ExperimentResults res;
    res.strategy_name = "Traditional";
    int prev_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    auto start_total = high_resolution_clock::now();

    auto start1 = high_resolution_clock::now();
    for (auto& net : cbn->l_local_networks) {
        if (net->local_scenes.empty()) {
            auto scenes = cbn->_generate_local_scenes(net);
            LocalNetwork::find_local_attractors_brute_force(net, scenes);
        }
        cbn->process_kind_signal(net);
    }
    cbn->_assign_global_indices_to_attractors();
    cbn->generate_attractor_dictionary();
    auto end1 = high_resolution_clock::now();
    res.p1_ms = duration<double, std::milli>(end1 - start1).count();

    auto start2 = high_resolution_clock::now();
    cbn->find_compatible_pairs();
    auto end2 = high_resolution_clock::now();
    res.p2_ms = duration<double, std::milli>(end2 - start2).count();

    auto start3 = high_resolution_clock::now();
    cbn->mount_stable_attractor_fields();
    auto end3 = high_resolution_clock::now();
    res.p3_ms = duration<double, std::milli>(end3 - start3).count();

    auto end_total = high_resolution_clock::now();
    res.total_ms = duration<double, std::milli>(end_total - start_total).count();
    res.max_rss_kb = get_max_rss();
    res.global_attractors_count = cbn->d_attractor_fields.size();
    res.success = true;
    omp_set_num_threads(prev_threads);
    return res;
}

// 2. SIMPLE PARALLEL EXPERIMENT
ExperimentResults SimpleParallelExperiment::run(std::shared_ptr<CBN> cbn) {
    ExperimentResults res;
    res.strategy_name = "SimpleParallel";
    auto start_total = high_resolution_clock::now();

    auto start1 = high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0; i < (int)cbn->l_local_networks.size(); ++i) {
        auto& net = cbn->l_local_networks.at(i);
        if (net->local_scenes.empty()) {
            auto scenes = cbn->_generate_local_scenes(net);
            LocalNetwork::find_local_attractors_brute_force(net, scenes);
        }
        cbn->process_kind_signal(net);
    }
    cbn->_assign_global_indices_to_attractors();
    cbn->generate_attractor_dictionary();
    auto end1 = high_resolution_clock::now();
    res.p1_ms = duration<double, std::milli>(end1 - start1).count();

    auto start2 = high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0; i < (int)cbn->l_directed_edges.size(); ++i) {
        auto& edge = cbn->l_directed_edges.at(i);
        for (int val : {0, 1}) {
            edge->d_comp_pairs_attractors_by_value[val].clear();
            auto dst_attractors = cbn->get_attractors_by_input_signal_value(edge->index_variable, val);
            auto& src_attractors = edge->d_out_value_to_attractor[val];
            edge->d_comp_pairs_attractors_by_value[val].reserve(src_attractors.size() * dst_attractors.size());
            for (auto& src_attr : src_attractors) {
                for (auto& dst_attr : dst_attractors) {
                    edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_attr->g_index});
                }
            }
        }
    }
    auto end2 = high_resolution_clock::now();
    res.p2_ms = duration<double, std::milli>(end2 - start2).count();

    auto start3 = high_resolution_clock::now();
    std::vector<std::vector<int>> fields;
    for (auto& edge : cbn->l_directed_edges) {
        for (int val : {0, 1}) {
            for (auto& pair : edge->d_comp_pairs_attractors_by_value[val]) {
                std::vector<int> f = {pair.first, pair.second};
                std::sort(f.begin(), f.end());
                fields.push_back(f);
            }
        }
    }

    if (!fields.empty()) {
        auto final_fields = merge_field_lists_sequentially(fields);
        cbn->d_attractor_fields.clear();
        int count = 1;
        for (auto& f : final_fields) cbn->d_attractor_fields[count++] = f;
    }
    auto end3 = high_resolution_clock::now();
    res.p3_ms = duration<double, std::milli>(end3 - start3).count();
    auto end_total = high_resolution_clock::now();
    res.total_ms = duration<double, std::milli>(end_total - start_total).count();
    res.max_rss_kb = get_max_rss();
    res.global_attractors_count = cbn->d_attractor_fields.size();
    res.success = true;
    return res;
}

// 3. ADVANCED PARALLEL EXPERIMENT
ExperimentResults AdvancedParallelExperiment::run(std::shared_ptr<CBN> cbn) {
    ExperimentResults res;
    res.strategy_name = "AdvancedParallel";
    auto start_total = high_resolution_clock::now();

    auto start1 = high_resolution_clock::now();
    struct WeightedNet { long weight; std::shared_ptr<LocalNetwork> net; };
    std::vector<WeightedNet> weighted_nets;
    for (auto& net : cbn->l_local_networks) {
        long w = (long)net->internal_variables.size() * (1L << net->input_signals.size());
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
        if (tid < (int)buckets.size()) {
            for (auto& net : buckets[tid]) {
                if (net->local_scenes.empty()) {
                    auto scenes = cbn->_generate_local_scenes(net);
                    LocalNetwork::find_local_attractors_brute_force(net, scenes);
                }
                cbn->process_kind_signal(net);
            }
        }
    }
    cbn->_assign_global_indices_to_attractors();
    cbn->generate_attractor_dictionary();
    auto end1 = high_resolution_clock::now();
    res.p1_ms = duration<double, std::milli>(end1 - start1).count();

    auto start2 = high_resolution_clock::now();
    std::map<int, std::vector<int>> var_to_gindices[2];
    for (auto& net : cbn->l_local_networks) {
        for (auto& scene : net->local_scenes) {
            for (size_t i = 0; i < scene->l_index_signals.size(); ++i) {
                int var_idx = scene->l_index_signals.at(i);
                int val = scene->l_values.at(0).at(i) - '0';
                if (val == 0 || val == 1) {
                    for (auto& attr : scene->l_attractors) var_to_gindices[val][var_idx].push_back(attr->g_index);
                }
            }
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < (int)cbn->l_directed_edges.size(); ++i) {
        auto& edge = cbn->l_directed_edges.at(i);
        for (int val : {0, 1}) {
            edge->d_comp_pairs_attractors_by_value[val].clear();
            auto& dst_g_indices = var_to_gindices[val][edge->index_variable];
            auto& src_attractors = edge->d_out_value_to_attractor[val];
            edge->d_comp_pairs_attractors_by_value[val].reserve(src_attractors.size() * dst_g_indices.size());
            for (auto& src_attr : src_attractors) {
                for (int dst_g_idx : dst_g_indices) edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_g_idx});
            }
        }
    }
    auto end2 = high_resolution_clock::now();
    res.p2_ms = duration<double, std::milli>(end2 - start2).count();

    auto start3 = high_resolution_clock::now();
    std::vector<std::vector<int>> initial_fields;
    for (auto& edge : cbn->l_directed_edges) {
        for (int val : {0, 1}) {
            for (auto& pair : edge->d_comp_pairs_attractors_by_value[val]) {
                std::vector<int> f = {pair.first, pair.second};
                std::sort(f.begin(), f.end());
                initial_fields.push_back(f);
            }
        }
    }

    std::function<std::vector<std::vector<int>>(std::vector<std::vector<int>>)> async_reduce;
    async_reduce = [&](std::vector<std::vector<int>> f_list) -> std::vector<std::vector<int>> {
        if (f_list.size() <= 100) return merge_field_lists_sequentially(f_list);
        size_t mid = f_list.size() / 2;
        std::vector<std::vector<int>> left_chunk(f_list.begin(), f_list.begin() + mid);
        std::vector<std::vector<int>> right_chunk(f_list.begin() + mid, f_list.end());

        auto left_f = std::async(std::launch::async, [&async_reduce, left_chunk]() {
            return async_reduce(left_chunk);
        });
        auto right_res = async_reduce(right_chunk);
        auto left_res = left_f.get();

        std::vector<std::vector<int>> combined = std::move(left_res);
        combined.insert(combined.end(), std::make_move_iterator(right_res.begin()), std::make_move_iterator(right_res.end()));
        return merge_field_lists_sequentially(std::move(combined));
    };

    if(!initial_fields.empty()) {
        auto final_fields = async_reduce(initial_fields);
        cbn->d_attractor_fields.clear();
        int count = 1;
        for (auto& f : final_fields) cbn->d_attractor_fields[count++] = f;
    }
    auto end3 = high_resolution_clock::now();
    res.p3_ms = duration<double, std::milli>(end3 - start3).count();
    auto end_total = high_resolution_clock::now();
    res.total_ms = duration<double, std::milli>(end_total - start_total).count();
    res.max_rss_kb = get_max_rss();
    res.global_attractors_count = cbn->d_attractor_fields.size();
    res.success = true;
    return res;
}

} // namespace cbnetwork
