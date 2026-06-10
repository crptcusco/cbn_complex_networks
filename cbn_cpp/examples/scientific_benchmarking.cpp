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

using namespace cbnetwork;
namespace fs = std::filesystem;

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
 * Persists sample data (fields and topology) to specific files.
 */
void persist_sample(int sample_id, std::shared_ptr<CBN> cbn, const ExperimentResults& res, const std::string& output_dir) {
    fs::create_directories(output_dir);

    // Persist Fields JSON: cbn_sample_{ID}_{STRATEGIA}_fields.json
    std::string fields_file = output_dir + "/cbn_sample_" + std::to_string(sample_id) + "_" + res.strategy_name + "_fields.json";
    std::ofstream f_out(fields_file);
    if (f_out.is_open()) {
        f_out << "{\n";
        f_out << "  \"sample_id\": " << sample_id << ",\n";
        f_out << "  \"strategy\": \"" << res.strategy_name << "\",\n";
        f_out << "  \"p1_ms\": " << res.p1_ms << ",\n";
        f_out << "  \"p2_ms\": " << res.p2_ms << ",\n";
        f_out << "  \"p3_ms\": " << res.p3_ms << ",\n";
        f_out << "  \"fields_count\": " << cbn->d_attractor_fields.size() << ",\n";
        f_out << "  \"fields\": [\n";
        size_t f_idx = 0;
        for (auto const& [id, field] : cbn->d_attractor_fields) {
            f_out << "    [";
            for (size_t i = 0; i < field.size(); ++i) {
                f_out << field[i] << (i == field.size() - 1 ? "" : ", ");
            }
            f_out << "]" << (++f_idx == cbn->d_attractor_fields.size() ? "" : ",\n");
        }
        f_out << "\n  ]\n";
        f_out << "}\n";
        f_out.close();
    } else {
        std::cerr << "[Error] Could not write to " << fields_file << std::endl;
    }

    // Persist Topology CSV: cbn_sample_{ID}_{STRATEGIA}_topology.csv
    std::string topo_file = output_dir + "/cbn_sample_" + std::to_string(sample_id) + "_" + res.strategy_name + "_topology.csv";
    std::ofstream t_out(topo_file);
    if (t_out.is_open()) {
        t_out << "source,target\n";
        for (auto& edge : cbn->l_directed_edges) {
            t_out << edge->output_local_network << "," << edge->input_local_network << "\n";
        }
        t_out.close();
    } else {
        std::cerr << "[Error] Could not write to " << topo_file << std::endl;
    }
}

int main(int argc, char** argv) {
    // Default values
    int n_samples = 3;
    int n_networks = 6;
    int n_vars = 10;
    int topology = 1;
    std::string output_file = "benchmark_results.csv";
    std::string output_dir = "phd_experiments/output_data";

    // CLI Parsing
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

    // Define strategies to test
    struct StrategyInfo {
        std::string name;
        std::unique_ptr<ExperimentStrategy> strat;
    };
    std::vector<StrategyInfo> strategies;
    strategies.push_back({"Traditional", std::make_unique<TraditionalExperiment>()});
    strategies.push_back({"SimpleParallel", std::make_unique<SimpleParallelExperiment>()});
    strategies.push_back({"AdvancedParallel", std::make_unique<AdvancedParallelExperiment>()});

    for (int s = 0; s < n_samples; ++s) {
        int sample_id = s + 1;
        std::cout << "\n[*] Processing Sample " << sample_id << "/" << n_samples << "..." << std::endl;

        for (auto& s_info : strategies) {
            try {
                // Generate a fresh CBN for each strategy to ensure fair comparison
                // and to avoid side-effects from previous runs (like cached attractors).
                auto cbn = CBN::cbn_generator(topology, n_networks, n_vars, 1, 1);

                if (!cbn) {
                    throw std::runtime_error("CBN Generator returned null for sample " + std::to_string(sample_id));
                }

                auto results = s_info.strat->run(cbn);

                // Logging and Persistence
                log_to_csv(output_file, results, sample_id, topology, n_networks, n_vars);
                persist_sample(sample_id, cbn, results, output_dir);

                std::cout << "  - " << std::left << std::setw(20) << results.strategy_name
                          << " | Fields: " << std::setw(5) << results.global_attractors_count
                          << " | Total: " << std::fixed << std::setprecision(2) << results.total_ms << " ms [OK]" << std::endl;

            } catch (const std::exception& e) {
                std::cerr << "  [!] Error in Sample " << sample_id << " (" << s_info.name << "): " << e.what() << std::endl;
                // We continue with the next strategy/sample
            } catch (...) {
                std::cerr << "  [!] Unknown error in Sample " << sample_id << " (" << s_info.name << ")" << std::endl;
            }
        }
    }

    std::cout << "\n==========================================================" << std::endl;
    std::cout << " Benchmarking completed. Results saved to " << output_file << std::endl;
    std::cout << "==========================================================" << std::endl;

    return 0;
}
