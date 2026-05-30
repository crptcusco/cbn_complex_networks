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

/**
 * Utility to track peak memory usage.
 */
static long get_max_rss() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

// ----------------------------------------------------------------------------
// 1. TRADITIONAL EXPERIMENT (Strict Baseline)
// ----------------------------------------------------------------------------
ExperimentResults TraditionalExperiment::run(std::shared_ptr<CBN> cbn) {
    ExperimentResults res;
    res.strategy_name = "Traditional";

    // Force strictly single-threaded execution for this trial
    int prev_threads = omp_get_max_threads();
    omp_set_num_threads(1);

    auto start_total = high_resolution_clock::now();

    // Phase 1: Sequential Attractors
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

    // Phase 2: Sequential Compatible Pairs
    auto start2 = high_resolution_clock::now();
    cbn->find_compatible_pairs();
    auto end2 = high_resolution_clock::now();
    res.p2_ms = duration<double, std::milli>(end2 - start2).count();

    // Phase 3: Sequential Attractor Fields
    auto start3 = high_resolution_clock::now();
    cbn->mount_stable_attractor_fields();
    auto end3 = high_resolution_clock::now();
    res.p3_ms = duration<double, std::milli>(end3 - start3).count();

    auto end_total = high_resolution_clock::now();
    res.total_ms = duration<double, std::milli>(end_total - start_total).count();
    res.max_rss_kb = get_max_rss();
    res.global_attractors_count = cbn->d_attractor_fields.size();
    res.success = true;

    // Restore previous thread count
    omp_set_num_threads(prev_threads);
    return res;
}

// ----------------------------------------------------------------------------
// 2. SIMPLE PARALLEL EXPERIMENT (Flat OpenMP)
// ----------------------------------------------------------------------------
ExperimentResults SimpleParallelExperiment::run(std::shared_ptr<CBN> cbn) {
    ExperimentResults res;
    res.strategy_name = "SimpleParallel";

    auto start_total = high_resolution_clock::now();

    // Phase 1: Flat Parallel Attractors
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

    // Phase 2: Flat Parallel Compatible Pairs
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

    // Phase 3: Flat Parallel Fields (Iterative merge with lock-free check)
    auto start3 = high_resolution_clock::now();
    std::vector<std::set<int>> fields;
    for (auto& edge : cbn->l_directed_edges) {
        for (int val : {0, 1}) {
            for (auto& pair : edge->d_comp_pairs_attractors_by_value[val]) {
                fields.push_back({pair.first, pair.second});
            }
        }
    }

    if (!fields.empty()) {
        bool global_merged = true;
        while (global_merged) {
            global_merged = false;
            for (size_t i = 0; i < fields.size(); ++i) {
                #pragma omp parallel for schedule(dynamic) shared(global_merged)
                for (int j = i + 1; j < (int)fields.size(); ++j) {
                    if (global_merged) continue;
                    bool share = false;
                    for (int attr : fields.at(i)) {
                        if (fields.at(j).count(attr)) { share = true; break; }
                    }
                    if (share) {
                        #pragma omp critical
                        {
                            if (!global_merged && i < fields.size() && j < (int)fields.size()) {
                                fields.at(i).insert(fields.at(j).begin(), fields.at(j).end());
                                fields.erase(fields.begin() + j);
                                global_merged = true;
                            }
                        }
                    }
                }
                if (global_merged) break;
            }
        }
        cbn->d_attractor_fields.clear();
        for (size_t i = 0; i < fields.size(); ++i) {
            cbn->d_attractor_fields[i + 1] = std::vector<int>(fields.at(i).begin(), fields.at(i).end());
        }
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

// ----------------------------------------------------------------------------
// 3. ADVANCED PARALLEL EXPERIMENT (Heuristic Balance & Asynchronous)
// ----------------------------------------------------------------------------
ExperimentResults AdvancedParallelExperiment::run(std::shared_ptr<CBN> cbn) {
    ExperimentResults res;
    res.strategy_name = "AdvancedParallel";

    auto start_total = high_resolution_clock::now();

    // Phase 1: Heuristic Bucket Balancing (V * 2^C)
    auto start1 = high_resolution_clock::now();
    struct WeightedNet {
        long weight;
        std::shared_ptr<LocalNetwork> net;
    };
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
        int best_idx = 0;
        long min_w = bucket_weights[0];
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

    // Phase 2: Variable-Indexed Map (O(1) equivalent)
    auto start2 = high_resolution_clock::now();
    std::map<int, std::vector<int>> var_to_gindices[2];
    for (auto& net : cbn->l_local_networks) {
        for (auto& scene : net->local_scenes) {
            for (size_t i = 0; i < scene->l_index_signals.size(); ++i) {
                int var_idx = scene->l_index_signals.at(i);
                int val = scene->l_values.at(0).at(i) - '0';
                if (val == 0 || val == 1) {
                    for (auto& attr : scene->l_attractors) {
                        var_to_gindices[val][var_idx].push_back(attr->g_index);
                    }
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
                for (int dst_g_idx : dst_g_indices) {
                    edge->d_comp_pairs_attractors_by_value[val].push_back({src_attr->g_index, dst_g_idx});
                }
            }
        }
    }
    auto end2 = high_resolution_clock::now();
    res.p2_ms = duration<double, std::milli>(end2 - start2).count();

    // Phase 3: Task-Based Parallelism (Binary Reduction Tree with std::async)
    auto start3 = high_resolution_clock::now();
    std::vector<std::set<int>> initial_fields;
    for (auto& edge : cbn->l_directed_edges) {
        for (int val : {0, 1}) {
            for (auto& pair : edge->d_comp_pairs_attractors_by_value[val]) {
                initial_fields.push_back({pair.first, pair.second});
            }
        }
    }

    std::function<std::vector<std::set<int>>(std::vector<std::set<int>>)> async_reduce;
    async_reduce = [&](std::vector<std::set<int>> f_list) -> std::vector<std::set<int>> {
        if (f_list.size() <= 20) {
            bool merged = true;
            while(merged) {
                merged = false;
                for(size_t i=0; i<f_list.size(); ++i) {
                    for(size_t j=i+1; j<f_list.size(); ++j) {
                        bool share = false;
                        for(int a : f_list[i]) if(f_list[j].count(a)) { share=true; break; }
                        if(share) { f_list[i].insert(f_list[j].begin(), f_list[j].end()); f_list.erase(f_list.begin()+j); merged=true; break; }
                    }
                    if(merged) break;
                }
            }
            return f_list;
        }
        size_t mid = f_list.size() / 2;
        auto left_f = std::async(std::launch::async, async_reduce, std::vector<std::set<int>>(f_list.begin(), f_list.begin()+mid));
        auto right_res = async_reduce(std::vector<std::set<int>>(f_list.begin()+mid, f_list.end()));
        auto left_res = left_f.get();
        for(auto& rf : right_res) {
            bool fm = false;
            for(auto& lf : left_res) {
                bool sh = false;
                for(int a : rf) if(lf.count(a)) { sh=true; break; }
                if(sh) { lf.insert(rf.begin(), rf.end()); fm=true; break; }
            }
            if(!fm) left_res.push_back(rf);
        }
        return left_res;
    };

    if(!initial_fields.empty()) {
        auto final_fields = async_reduce(initial_fields);
        cbn->d_attractor_fields.clear();
        for (size_t i = 0; i < final_fields.size(); ++i) {
            cbn->d_attractor_fields[i + 1] = std::vector<int>(final_fields.at(i).begin(), final_fields.at(i).end());
        }
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
