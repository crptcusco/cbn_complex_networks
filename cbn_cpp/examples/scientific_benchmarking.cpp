#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <iomanip>
#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/experiment_strategies.hpp"

using namespace cbnetwork;

void log_to_csv(const std::string& filename, const ExperimentResults& res, int sample_id, int topology, int networks, int vars) {
    std::ofstream out(filename, std::ios_base::app);
    out << sample_id << "," << res.strategy_name << "," << topology << "," << networks << "," << vars << ","
        << res.p1_ms << "," << res.p2_ms << "," << res.p3_ms << ","
        << res.max_rss_kb << "," << res.global_attractors_count << "\n";
    out.flush();
    out.close();
}

int main(int argc, char** argv) {
    int n_samples = 5;
    int n_networks = 5;
    int n_vars = 6;
    int topology = 1;
    std::string output_file = "benchmark_results.csv";

    // Simple arg parsing
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

    std::vector<std::unique_ptr<ExperimentStrategy>> strategies;
    strategies.push_back(std::make_unique<TraditionalExperiment>());
    strategies.push_back(std::make_unique<SimpleParallelExperiment>());
    strategies.push_back(std::make_unique<AdvancedParallelExperiment>());

    for (int s = 0; s < n_samples; ++s) {
        std::cout << "Starting Sample " << s + 1 << "/" << n_samples << "..." << std::endl;

        // Generate once
        auto cbn_master = CBN::cbn_generator(topology, n_networks, n_vars, 2, 2);

        for (auto& strat : strategies) {
            // Since Step 1/2 modify the object state (attractors/pairs),
            // we need to either clone or re-generate.
            // Deep cloning is not implemented. Re-generating with same master topology but new dynamics.
            // For rigorous benchmarking, they should see identical logic.
            // Let's pass the master to Traditional, then reset for others or just re-generate.
            // Actually, generating a new one with same parameters is common in synthetic benchmarks.

            auto cbn = CBN::cbn_generator(topology, n_networks, n_vars, 2, 2);
            auto results = strat->run(cbn);
            log_to_csv(output_file, results, s + 1, topology, n_networks, n_vars);
            std::cout << "  " << std::left << std::setw(20) << results.strategy_name
                      << " completed in " << std::fixed << std::setprecision(2) << results.total_ms << " ms" << std::endl;
        }
    }

    return 0;
}
