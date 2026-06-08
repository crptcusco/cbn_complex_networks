#include "cbnetwork/campaign_manager.hpp"
#include "cbnetwork/network_factory.hpp"
#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <sys/resource.h>

namespace cbnetwork {

namespace fs = std::filesystem;
using json = nlohmann::json;

std::vector<ExperimentConfig> CampaignManager::load_campaign_config(const std::string& filepath) {
    std::vector<ExperimentConfig> experiments;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Could not open campaign config: " << filepath << std::endl;
        return experiments;
    }

    json config;
    file >> config;

    for (const auto& item : config["experiments"]) {
        ExperimentConfig exp;
        exp.experiment_id = item["id"];
        exp.topology_type = item["topology"];
        exp.num_networks = item["num_networks"];
        exp.vars_per_network = item["vars_per_network"];
        exp.input_variables = item["input_variables"];
        exp.output_variables = item["output_variables"];
        exp.max_clauses = item.value("max_clauses", 3);
        exp.max_literals = item.value("max_literals", 2);
        exp.coupling_strategy_name = item["coupling_strategy"];
        exp.threshold_value = item.value("threshold", 0);
        experiments.push_back(exp);
    }

    return experiments;
}

void CampaignManager::run_campaign(const std::vector<ExperimentConfig>& experiments) {
    for (const auto& config : experiments) {
        std::cout << "Starting Experiment: " << config.experiment_id << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        auto cbn = NetworkFactory::create_cbn(config);

        cbn->find_local_attractors_parallel();
        cbn->find_compatible_pairs_parallel();
        cbn->mount_stable_attractor_fields_parallel();

        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration<double>(end - start).count();

        export_experiment_data(config, cbn, duration);
        std::cout << "Completed Experiment: " << config.experiment_id << " in " << duration << "s" << std::endl;
    }
}

void CampaignManager::export_experiment_data(const ExperimentConfig& config,
                                           std::shared_ptr<CBN> cbn,
                                           double execution_time_s) {
    std::string dir_path = "experiments/" + config.experiment_id;
    fs::create_directories(dir_path);

    // 1. Network File
    cbn->save_network_to_json(dir_path + "/" + config.experiment_id + "_network.json");

    // 2. Dynamics File
    json j_dynamics;

    json j_local = json::array();
    for (const auto& net : cbn->l_local_networks) {
        json j_net;
        j_net["network_index"] = net->index;
        j_net["attractor_count"] = net->attractor_count;
        j_local.push_back(j_net);
    }
    j_dynamics["local_attractors"] = j_local;

    json j_edges = json::array();
    for (const auto& edge : cbn->l_directed_edges) {
        json j_edge;
        j_edge["edge_index"] = edge->index;
        j_edge["compat_pairs_0"] = edge->d_comp_pairs_attractors_by_value[0].size();
        j_edge["compat_pairs_1"] = edge->d_comp_pairs_attractors_by_value[1].size();
        j_edges.push_back(j_edge);
    }
    j_dynamics["compatible_pairs"] = j_edges;

    json j_fields = json::array();
    for (const auto& [id, field] : cbn->d_attractor_fields) {
        json j_f;
        j_f["field_id"] = id;
        j_f["attractors"] = field;
        j_fields.push_back(j_f);
    }
    j_dynamics["global_attractor_fields"] = j_fields;

    std::ofstream dyn_file(dir_path + "/" + config.experiment_id + "_dynamics.json");
    dyn_file << j_dynamics.dump(4) << std::endl;

    // 3. Performance File
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long memory_kb = usage.ru_maxrss; // On Linux this is in KB

    json j_perf;
    j_perf["execution_time_s"] = execution_time_s;
    j_perf["memory_usage_kb"] = memory_kb;
    j_perf["canalization_density"] = cbn->calculate_canalization_density();
    j_perf["n_local_attractors"] = cbn->get_n_local_attractors();
    j_perf["n_attractor_fields"] = cbn->get_n_attractor_fields();

    std::ofstream perf_file(dir_path + "/" + config.experiment_id + "_performance.json");
    perf_file << j_perf.dump(4) << std::endl;
}

} // namespace cbnetwork
