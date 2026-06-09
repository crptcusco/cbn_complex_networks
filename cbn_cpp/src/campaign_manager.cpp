#include "cbnetwork/campaign_manager.hpp"
#include "cbnetwork/network_factory.hpp"
#include "cbnetwork/experiment_strategies.hpp"
#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath>

namespace cbnetwork {

namespace fs = std::filesystem;
using json = nlohmann::json;

void CampaignManager::run_campaign(const std::string& config_path) {
    std::ifstream file(config_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open campaign config file: " << config_path << std::endl;
        return;
    }

    json j_config;
    file >> j_config;

    if (!j_config.contains("experiments") || !j_config["experiments"].is_array()) {
        std::cerr << "Error: Invalid campaign config format. Missing 'experiments' array." << std::endl;
        return;
    }

    fs::create_directories("experimentos");

    for (const auto& j_exp : j_config["experiments"]) {
        ExperimentConfig config;
        config.id = j_exp["id"].get<std::string>();
        config.topology_type = j_exp["topology_type"].get<int>();
        config.n_local_networks = j_exp["n_local_networks"].get<int>();
        config.n_vars_network = j_exp["n_vars_network"].get<int>();
        config.n_input_variables = j_exp["n_input_variables"].get<int>();
        config.n_output_variables = j_exp["n_output_variables"].get<int>();
        config.n_max_of_clauses = j_exp.value("n_max_of_clauses", 2);
        config.n_max_of_literals = j_exp.value("n_max_of_literals", 3);
        config.n_edges = j_exp.value("n_edges", -1);
        config.coupling_type = j_exp["coupling_type"].get<std::string>();
        config.coupling_param = j_exp.value("coupling_param", 0);

        std::cout << "Executing Experiment: " << config.id << " ..." << std::endl;

        std::string exp_dir = "experimentos/" + config.id;
        fs::create_directories(exp_dir);

        auto cbn = NetworkFactory::create_from_config(config);

        AdvancedParallelExperiment strategy;
        ExperimentResults results = strategy.run(cbn);

        if (results.success) {
            export_experiment_data(exp_dir + "/" + config.id, config, cbn, results);
            std::cout << "Experiment " << config.id << " completed successfully." << std::endl;
        } else {
            std::cerr << "Experiment " << config.id << " failed." << std::endl;
        }
    }
}

void CampaignManager::export_experiment_data(const std::string& base_path,
                                            const ExperimentConfig& config,
                                            std::shared_ptr<CBN> cbn,
                                            const ExperimentResults& results) {
    // 1. [id]_network.json
    cbn->save_network_to_json(base_path + "_network.json");

    // 2. [id]_dynamics.json
    json j_dynamics;

    // Compatible Pairs (from edges)
    json j_edges_pairs = json::array();
    for (const auto& edge : cbn->l_directed_edges) {
        json j_edge;
        j_edge["edge_index"] = edge->index;
        j_edge["signal_var"] = edge->index_variable;
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

    // HIARARCHICAL STRUCTURE: global_scenes
    json j_global_scenes = json::array();

    int num_global_edges = cbn->l_directed_edges.size();
    if (num_global_edges > 0 && num_global_edges <= 20) {
        long long total_scenes = 1LL << num_global_edges;

        for (long long i = 0; i < total_scenes; ++i) {
            json j_scene;
            j_scene["global_scene_index"] = i + 1;

            std::vector<int> edge_states;
            std::map<int, int> edge_val_map; // edge_index -> state
            for (int j = 0; j < num_global_edges; ++j) {
                int state = (i >> j) & 1;
                edge_states.push_back(state);
                edge_val_map[cbn->l_directed_edges[j]->index] = state;
            }
            j_scene["global_edges_state"] = edge_states;

            json j_local_attractors_by_module = json::array();

            for (const auto& net : cbn->l_local_networks) {
                json j_module;
                j_module["module_id"] = net->index;
                json j_attrs = json::array();

                // Find attractors compatible with this global scene
                for (const auto& scene : net->local_scenes) {
                    bool match = true;
                    for (size_t k = 0; k < scene->l_index_signals.size(); ++k) {
                        int signal_var = scene->l_index_signals[k];
                        int expected_val = scene->l_values[0][k] - '0';

                        // Find which edge corresponds to this signal_var
                        bool found_edge = false;
                        for (const auto& edge : cbn->l_directed_edges) {
                            if (edge->index_variable == signal_var) {
                                if (edge_val_map[edge->index] != expected_val) {
                                    match = false;
                                }
                                found_edge = true;
                                break;
                            }
                        }
                        if (!match) break;
                    }

                    if (match) {
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
                }
                j_module["attractors"] = j_attrs;
                j_local_attractors_by_module.push_back(j_module);
            }
            j_scene["local_attractors_by_module"] = j_local_attractors_by_module;
            j_global_scenes.push_back(j_scene);
        }
    }
    j_dynamics["global_scenes"] = j_global_scenes;

    // Global Attractor Fields
    json j_fields = json::array();
    for (auto const& [id, field] : cbn->d_attractor_fields) {
        json j_field;
        j_field["field_id"] = id;
        j_field["attractor_indices"] = field;
        j_fields.push_back(j_field);
    }
    j_dynamics["global_attractor_fields"] = j_fields;

    std::ofstream dynamics_out(base_path + "_dynamics.json");
    dynamics_out << j_dynamics.dump(4) << std::endl;

    // 3. [id]_performance.json
    json j_performance;
    j_performance["execution_time_ms"] = {
        {"total", results.total_ms},
        {"phase1_local", results.p1_ms},
        {"phase2_pairs", results.p2_ms},
        {"phase3_fields", results.p3_ms}
    };
    j_performance["memory_usage_kb"] = results.max_rss_kb;
    j_performance["canalization_density"] = cbn->calculate_canalization_density();
    j_performance["counts"] = {
        {"local_attractors", cbn->get_n_local_attractors()},
        {"global_fields", results.global_attractors_count}
    };

    std::ofstream performance_out(base_path + "_performance.json");
    performance_out << j_performance.dump(4) << std::endl;
}

} // namespace cbnetwork
