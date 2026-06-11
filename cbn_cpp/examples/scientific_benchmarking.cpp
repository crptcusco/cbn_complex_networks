#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/experiment_strategies.hpp"
#include "nlohmann/json.hpp"

using namespace cbnetwork;
namespace fs = std::filesystem;
using json = nlohmann::json;

/**
 * Logs metrics to a global CSV file.
 */
void log_to_csv(const std::string& filename, const ExperimentResults& res, int sample_id, int topology, int networks, int vars) {
    bool file_exists = fs::exists(filename);
    std::ofstream out(filename, std::ios_base::app);
    if (!out.is_open()) {
        std::cerr << "[Error] Could not open log file: " << filename << std::endl;
        return;
    }

    if (!file_exists) {
        out << "sample_id,strategy,topology,n_networks,n_vars,p1_ms,p2_ms,p3_ms,max_rss_kb,fields\n";
    }

    out << sample_id << "," << res.strategy_name << "," << topology << "," << networks << "," << vars << ","
        << std::fixed << std::setprecision(4) << res.p1_ms << ","
        << std::fixed << std::setprecision(4) << res.p2_ms << ","
        << std::fixed << std::setprecision(4) << res.p3_ms << ","
        << res.max_rss_kb << "," << res.global_attractors_count << "\n";

    out.flush();
    out.close();
}

/**
 * Exports network structure to JSON.
 */
void export_network_structure_json(int sample_id, std::shared_ptr<CBN> cbn, const std::string& output_dir) {
    std::string filepath = output_dir + "/cbn_sample_" + std::to_string(sample_id) + "_network_structure.json";
    cbn->save_network_to_json(filepath);
}

/**
 * Exports dynamics data to JSON.
 */
void export_dynamics_json(int sample_id, const std::string& strategy_name, std::shared_ptr<CBN> cbn, const ExperimentResults& res, const std::string& output_dir) {
    std::string filepath = output_dir + "/cbn_sample_" + std::to_string(sample_id) + "_" + strategy_name + "_dynamics.json";

    json j_dynamics;
    j_dynamics["sample_id"] = sample_id;
    j_dynamics["strategy"] = strategy_name;
    j_dynamics["performance"] = {
        {"p1_ms", res.p1_ms},
        {"p2_ms", res.p2_ms},
        {"p3_ms", res.p3_ms},
        {"total_ms", res.total_ms}
    };

    // Local Attractors
    json j_local_attractors = json::array();
    for (const auto& net : cbn->l_local_networks) {
        json j_net_attractors;
        j_net_attractors["network_index"] = net->index;
        json j_attrs = json::array();
        for (const auto& scene : net->local_scenes) {
            for (const auto& attr : scene->l_attractors) {
                json j_attr;
                j_attr["global_index"] = attr->g_index;
                j_attr["local_index"] = attr->l_index;
                j_attr["scene_index"] = scene->index;
                json j_states = json::array();
                for (const auto& state : attr->l_states) {
                    j_states.push_back(state->l_variable_values);
                }
                j_attr["states"] = j_states;
                j_attrs.push_back(j_attr);
            }
        }
        j_net_attractors["attractors"] = j_attrs;
        j_local_attractors.push_back(j_net_attractors);
    }
    j_dynamics["local_attractors"] = j_local_attractors;

    // Compatible Pairs
    json j_edges_pairs = json::array();
    for (const auto& edge : cbn->l_directed_edges) {
        json j_edge;
        j_edge["edge_index"] = edge->index;
        json j_vals = json::array();
        for (int val : {0, 1}) {
            json j_val_pairs;
            j_val_pairs["signal_value"] = val;
            j_val_pairs["pairs"] = edge->d_comp_pairs_attractors_by_value[val];
            j_vals.push_back(j_val_pairs);
        }
        j_edge["compatible_pairs"] = j_vals;
        j_edges_pairs.push_back(j_edge);
    }
    j_dynamics["compatible_pairs_by_edge"] = j_edges_pairs;

    // Global Attractor Fields
    json j_fields = json::array();
    for (auto const& [id, field] : cbn->d_attractor_fields) {
        json j_field;
        j_field["field_id"] = id;
        j_field["attractor_indices"] = field;
        j_fields.push_back(j_field);
    }
    j_dynamics["global_attractor_fields"] = j_fields;

    std::ofstream out(filepath);
    if (out.is_open()) {
        out << j_dynamics.dump(4) << std::endl;
        out.close();
    } else {
        std::cerr << "[Error] Could not write to " << filepath << std::endl;
    }
}

/**
 * Exports topology to CSV.
 */
void export_topology_csv(int sample_id, const std::string& strategy_name, std::shared_ptr<CBN> cbn, const std::string& output_dir) {
    std::string filepath = output_dir + "/cbn_sample_" + std::to_string(sample_id) + "_" + strategy_name + "_topology.csv";
    std::ofstream out(filepath);
    if (out.is_open()) {
        out << "source,target\n";
        for (auto& edge : cbn->l_directed_edges) {
            out << edge->output_local_network << "," << edge->input_local_network << "\n";
        }
        out.close();
    } else {
        std::cerr << "[Error] Could not write to " << filepath << std::endl;
    }
}

int main(int argc, char** argv) {
    int n_samples = 3;
    int n_networks = 6;
    int n_vars = 10;
    int topology = 1;
    std::string output_file = "benchmark_results.csv";
    std::string output_dir = "phd_experiments/output_data";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        try {
            if (arg == "--samples" && i + 1 < argc) n_samples = std::stoi(argv[++i]);
            else if ((arg == "--networks" || arg == "--nodes") && i + 1 < argc) n_networks = std::stoi(argv[++i]);
            else if (arg == "--vars" && i + 1 < argc) n_vars = std::stoi(argv[++i]);
            else if (arg == "--topology" && i + 1 < argc) topology = std::stoi(argv[++i]);
            else if (arg == "--output" && i + 1 < argc) output_file = argv[++i];
            else if (arg == "--dir" && i + 1 < argc) output_dir = argv[++i];
        } catch (const std::exception& e) {
            std::cerr << "[Error] Parsing argument " << arg << ": " << e.what() << std::endl;
            return 1;
        }
    }

    fs::create_directories(output_dir);

    std::cout << "==========================================================" << std::endl;
    std::cout << " SCIENTIFIC BENCHMARKING PIPELINE " << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << " Samples:   " << n_samples << std::endl;
    std::cout << " Networks:  " << n_networks << std::endl;
    std::cout << " Variables: " << n_vars << std::endl;
    std::cout << " Topology:  " << topology << std::endl;
    std::cout << " Output:    " << output_file << std::endl;
    std::cout << " Dir:       " << output_dir << std::endl;
    std::cout << "==========================================================" << std::endl;

    std::vector<std::pair<std::string, std::unique_ptr<ExperimentStrategy>>> strategies;
    strategies.push_back({"Traditional", std::make_unique<TraditionalExperiment>()});
    strategies.push_back({"SimpleParallel", std::make_unique<SimpleParallelExperiment>()});
    strategies.push_back({"AdvancedParallel", std::make_unique<AdvancedParallelExperiment>()});

    for (int s = 0; s < n_samples; ++s) {
        int sample_id = s + 1;
        std::cout << "\n[*] Processing Sample " << sample_id << "/" << n_samples << "..." << std::endl;

        try {
            // Instantiate reference CBN per sample
            auto reference_cbn = CBN::cbn_generator(topology, n_networks, n_vars, 1, 1);
            if (!reference_cbn) throw std::runtime_error("CBN Generator failed.");

            // Export Structure (once per sample)
            export_network_structure_json(sample_id, reference_cbn, output_dir);

            for (auto& s_pair : strategies) {
                const std::string& s_name = s_pair.first;
                auto& s_strat = s_pair.second;

                try {
                    // DEEP COPY for absolute isolation and scientific parity
                    auto cbn = reference_cbn->clone();

                    auto results = s_strat->run(cbn);

                    // Logging and Persistence
                    log_to_csv(output_file, results, sample_id, topology, n_networks, n_vars);
                    export_dynamics_json(sample_id, s_name, cbn, results, output_dir);
                    export_topology_csv(sample_id, s_name, cbn, output_dir);

                    std::cout << "  - " << std::left << std::setw(20) << s_name
                              << " | Fields: " << std::setw(5) << results.global_attractors_count
                              << " | Total: " << std::fixed << std::setprecision(2) << results.total_ms << " ms [OK]" << std::endl;

                } catch (const std::exception& e) {
                    std::cerr << "  [!] Error in " << s_name << ": " << e.what() << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << " [!] Error generating Sample " << sample_id << ": " << e.what() << std::endl;
        }
    }

    std::cout << "\n==========================================================" << std::endl;
    std::cout << " Benchmarking completed. " << std::endl;
    std::cout << "==========================================================" << std::endl;

    return 0;
}
