#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <filesystem>
#include "nlohmann/json.hpp"
#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/localnetwork.hpp"
#include "cbnetwork/directededge.hpp"
#include "cbnetwork/internalvariable.hpp"
#include "cbnetwork/path_utils.hpp"

using json = nlohmann::json;
using namespace cbnetwork;
using namespace std::chrono;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <network_config.json> [--topology <id>] [--networks <n>] [--vars <n>] [--samples <n>]" << std::endl;
        return 1;
    }

    std::string config_file = argv[1];
    int topology = 0;
    int n_networks = 0;
    int n_vars = 0;
    int n_samples = 1;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--topology" && i + 1 < argc) topology = std::stoi(argv[++i]);
        else if (arg == "--networks" && i + 1 < argc) n_networks = std::stoi(argv[++i]);
        else if (arg == "--vars" && i + 1 < argc) n_vars = std::stoi(argv[++i]);
        else if (arg == "--samples" && i + 1 < argc) n_samples = std::stoi(argv[++i]);
    }

    std::ifstream i(config_file);
    if (!i.is_open()) {
        std::cerr << "Could not open file: " << config_file << std::endl;
        return 1;
    }

    json j;
    try {
        i >> j;
    } catch (json::parse_error& e) {
        std::cerr << "Parse error: " << e.what() << std::endl;
        return 1;
    }

    std::vector<std::shared_ptr<LocalNetwork>> networks;
    std::map<int, std::shared_ptr<LocalNetwork>> net_map;

    // Load Local Networks
    if (!j.contains("local_networks")) {
        std::cerr << "JSON missing 'local_networks'" << std::endl;
        return 1;
    }

    for (auto& net_j : j["local_networks"]) {
        int idx = net_j["index"];
        std::vector<int> internal_vars = net_j["internal_variables"].get<std::vector<int>>();
        auto net = std::make_shared<LocalNetwork>(idx, internal_vars);

        std::string logic_key = net_j.contains("descriptive_function_variables") ? "descriptive_function_variables" : "logic";
        for (auto& var_j : net_j[logic_key]) {
            int v_idx = var_j["index"];
            std::vector<std::vector<int>> cnf = var_j["cnf"].get<std::vector<std::vector<int>>>();
            auto var = std::make_shared<InternalVariable>(v_idx, cnf);
            net->descriptive_function_variables.push_back(var);
        }
        networks.push_back(net);
        net_map[idx] = net;
    }

    // Update info if not provided but in JSON (best effort)
    if (n_networks == 0) n_networks = networks.size();
    if (n_vars == 0 && !networks.empty()) n_vars = networks[0]->internal_variables.size();

    std::vector<std::shared_ptr<DirectedEdge>> edges;
    // Load Directed Edges
    if (j.contains("directed_edges")) {
        for (auto& edge_j : j["directed_edges"]) {
            int idx = edge_j["index"];

            // Handle both naming conventions
            int idx_var = edge_j.contains("index_variable") ? edge_j["index_variable"].get<int>() : edge_j["signal_var"].get<int>();
            int in_net = edge_j.contains("input_local_network") ? edge_j["input_local_network"].get<int>() : edge_j["target"].get<int>();
            int out_net = edge_j.contains("output_local_network") ? edge_j["output_local_network"].get<int>() : edge_j["source"].get<int>();
            std::vector<int> out_vars = edge_j.contains("output_variables") ? edge_j["output_variables"].get<std::vector<int>>() : edge_j["output_vars"].get<std::vector<int>>();
            std::string coup_func = edge_j.contains("coupling_function") ? edge_j["coupling_function"].get<std::string>() : edge_j["function"].get<std::string>();

            auto edge = std::make_shared<DirectedEdge>(idx, idx_var, in_net, out_net, out_vars, coup_func);

            if (edge_j.contains("true_table")) {
                std::map<std::string, std::string> tt = edge_j["true_table"].get<std::map<std::string, std::string>>();
                edge->true_table = tt;
            }

            edges.push_back(edge);
        }
    }

    // Connect input signals to networks
    for (auto& net : networks) {
        std::vector<std::shared_ptr<DirectedEdge>> inputs;
        for (auto& edge : edges) {
            if (edge->input_local_network == net->index) {
                inputs.push_back(edge);
            }
        }
        net->process_input_signals(inputs);
    }

    auto start_total = high_resolution_clock::now();

    CBN cbn(networks, edges);
    cbn.process_output_signals();

    // Step 1
    auto s1_start = high_resolution_clock::now();
    cbn.find_local_attractors_sequential();
    auto s1_end = high_resolution_clock::now();
    double s1_time = duration<double>(s1_end - s1_start).count();

    // Step 2
    auto s2_start = high_resolution_clock::now();
    cbn.find_compatible_pairs();
    auto s2_end = high_resolution_clock::now();
    double s2_time = duration<double>(s2_end - s2_start).count();

    // Step 3
    auto s3_start = high_resolution_clock::now();
    cbn.mount_stable_attractor_fields();
    auto s3_end = high_resolution_clock::now();
    double s3_time = duration<double>(s3_end - s3_start).count();

    auto end_total = high_resolution_clock::now();
    double total_time = duration<double>(end_total - start_total).count();

    // Collect metrics
    json res;
    res["status"] = "success";
    res["metrics"] = {
        {"step1_time", s1_time},
        {"step2_time", s2_time},
        {"step3_time", s3_time},
        {"total_time", total_time},
        {"local_attractors_count", 0},
        {"compatible_pairs_count", 0},
        {"attractor_fields_count", (int)cbn.d_attractor_fields.size()}
    };

    int total_local_attractors = 0;
    for(auto& net : cbn.l_local_networks) {
        total_local_attractors += net->attractor_count;
    }
    res["metrics"]["local_attractors_count"] = total_local_attractors;

    int total_pairs = 0;
    for(auto& edge : cbn.l_directed_edges) {
        total_pairs += edge->d_comp_pairs_attractors_by_value[0].size();
        total_pairs += edge->d_comp_pairs_attractors_by_value[1].size();
    }
    res["metrics"]["compatible_pairs_count"] = total_pairs;

    // Handle dynamic output directory
    if (topology > 0 || n_networks > 0 || n_vars > 0) {
        std::string output_dir = create_output_directory(topology, n_networks, n_vars, n_samples);
        std::string base_name = fs::path(config_file).stem().string();

        std::string metrics_file = output_dir + "/" + base_name + "_metrics.json";
        std::ofstream o_metrics(metrics_file);
        o_metrics << res.dump(4) << std::endl;

        std::string struct_file = output_dir + "/" + base_name + "_structure.json";
        cbn.save_network_to_json(struct_file);

        std::cout << "Results saved to: " << output_dir << std::endl;
    }

    std::cout << res.dump(4) << std::endl;

    return 0;
}
