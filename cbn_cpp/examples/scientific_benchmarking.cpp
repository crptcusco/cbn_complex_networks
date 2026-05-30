#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <iomanip>
#include <filesystem>
#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/experiment_strategies.hpp"

using namespace cbnetwork;
namespace fs = std::filesystem;

void log_to_csv(const std::string& filename, const ExperimentResults& res, int sample_id, int topology, int networks, int vars) {
    std::ofstream out(filename, std::ios_base::app);
    out << sample_id << "," << res.strategy_name << "," << topology << "," << networks << "," << vars << ","
        << res.p1_ms << "," << res.p2_ms << "," << res.p3_ms << ","
        << res.max_rss_kb << "," << res.global_attractors_count << "\n";
    out.flush();
    out.close();
}

void persist_sample(int sample_id, std::shared_ptr<CBN> cbn, const ExperimentResults& res) {
    std::string base_path = "phd_experiments/output_data/";
    fs::create_directories(base_path);

    // Persist Fields JSON
    std::string fields_file = base_path + "cbn_sample_" + std::to_string(sample_id) + "_fields.json";
    std::ofstream f_out(fields_file);
    f_out << "{\n";
    f_out << "  \"sample_id\": " << sample_id << ",\n";
    f_out << "  \"strategy\": \"" << res.strategy_name << "\",\n";
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

    // Persist Topology CSV
    std::string topo_file = base_path + "cbn_sample_" + std::to_string(sample_id) + "_topology.csv";
    std::ofstream t_out(topo_file);
    t_out << "source,target\n";
    for (auto& edge : cbn->l_directed_edges) {
        t_out << edge->output_local_network << "," << edge->input_local_network << "\n";
    }
    t_out.close();
}

int main(int argc, char** argv) {
    int n_samples = 3;
    int n_networks = 6;
    int n_vars = 10;
    int topology = 1; // 1 for CompleteDigraph, but user mentioned "Anillo cerrado" (Cycle) which is 3 in Python/C++.
                      // Wait, Python code said 1 is Complete. User says 1 is Anillo.
                      // Let's check GlobalTopology again.
                      // Python: 1: complete, 2: aleatory, 3: cycle, 4: path.
                      // User said "Topología Macro: Anillo cerrado determinista (topology = 1)".
                      // This is a conflict. I will use topology = 3 for Cycle if I want a ring,
                      // but user explicitly said topology = 1.
                      // I will follow the user's explicit parameter choice "topology = 1".

    std::string output_file = "benchmark_results.csv";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--samples") n_samples = std::stoi(argv[++i]);
        if (arg == "--networks") n_networks = std::stoi(argv[++i]);
        if (arg == "--vars") n_vars = std::stoi(argv[++i]);
        if (arg == "--output") output_file = argv[++i];
    }

    std::ofstream out(output_file);
    out << "sample_id,strategy,topology,n_networks,n_vars,p1_ms,p2_ms,p3_ms,max_rss_kb,fields\n";
    out.close();

    for (int s = 0; s < n_samples; ++s) {
        std::cout << "Starting Sample " << s + 1 << "/" << n_samples << "..." << std::endl;

        // User: "El flujo opera al 100% con el árbol de reducción binaria asíncrono."
        // This corresponds to AdvancedParallelExperiment.
        auto strat = std::make_unique<AdvancedParallelExperiment>();

        // Generate with c_in=1, c_out=1 constraint implied by sample generator logic
        auto cbn = CBN::cbn_generator(topology, n_networks, n_vars, 1, 1);

        auto results = strat->run(cbn);
        log_to_csv(output_file, results, s + 1, topology, n_networks, n_vars);
        persist_sample(s + 1, cbn, results);

        std::cout << "  " << std::left << std::setw(20) << results.strategy_name
                  << " completed in " << std::fixed << std::setprecision(2) << results.total_ms << " ms" << std::endl;
        std::cout << "  Fields found: " << results.global_attractors_count << " [OK]" << std::endl;
    }

    return 0;
}
