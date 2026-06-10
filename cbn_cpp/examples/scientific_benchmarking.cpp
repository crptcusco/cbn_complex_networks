#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/experiment_strategies.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace cbnetwork;
namespace fs = std::filesystem;

void log_to_csv(const std::string &filename, const ExperimentResults &res,
                int sample_id, int topology, int networks, int vars) {
  std::ofstream out(filename, std::ios_base::app);

  // EL DETECTOR DE FALLOS SILENCIOSOS
  if (!out.is_open()) {
    std::cerr << "ERROR FATAL: C++ no tiene permisos para escribir en "
              << filename << ". ¿El archivo está abierto en otro programa?\n";
    return;
  }

  out << sample_id << "," << res.strategy_name << "," << topology << ","
      << networks << "," << vars << "," << res.p1_ms << "," << res.p2_ms << ","
      << res.p3_ms << "," << res.max_rss_kb << ","
      << res.global_attractors_count << "\n";
  out.close();
}

void persist_sample(int sample_id, std::shared_ptr<CBN> cbn,
                    const ExperimentResults &res) {
  std::string base_path = "phd_experiments/output_data/";
  fs::create_directories(base_path);

  // Persist Fields JSON
  std::string fields_file =
      base_path + "cbn_sample_" + std::to_string(sample_id) + "_fields.json";
  std::ofstream f_out(fields_file);
  f_out << "{\n";
  f_out << "  \"sample_id\": " << sample_id << ",\n";
  f_out << "  \"strategy\": \"" << res.strategy_name << "\",\n";
  f_out << "  \"p3_ms\": " << res.p3_ms << ",\n";
  f_out << "  \"fields_count\": " << cbn->d_attractor_fields.size() << ",\n";
  f_out << "  \"fields\": [\n";
  size_t f_idx = 0;
  for (auto const &[id, field] : cbn->d_attractor_fields) {
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
  std::string topo_file =
      base_path + "cbn_sample_" + std::to_string(sample_id) + "_topology.csv";
  std::ofstream t_out(topo_file);
  t_out << "source,target\n";
  for (auto &edge : cbn->l_directed_edges) {
    t_out << edge->output_local_network << "," << edge->input_local_network
          << "\n";
  }
  t_out.close();
}

int main(int argc, char **argv) {
  int n_samples = 100;
  int n_networks = 11;
  int n_vars = 10;
  int topology = 4;
  std::string output_file = "benchmark_results.csv";

  // LECTURA DE ARGUMENTOS CORREGIDA
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--samples")
      n_samples = std::stoi(argv[++i]);
    else if (arg == "--networks" || arg == "--nodes")
      n_networks = std::stoi(argv[++i]);
    else if (arg == "--vars")
      n_vars = std::stoi(argv[++i]);
    else if (arg == "--output")
      output_file = argv[++i];
    else if (arg == "--topology")
      topology = std::stoi(argv[++i]);
  }

  // Abrimos el archivo una sola vez al inicio para crear/limpiar el CSV
  std::ofstream out(output_file);
  out << "sample_id,strategy,topology,n_networks,n_vars,p1_ms,p2_ms,p3_ms,max_"
         "rss_kb,fields\n";
  out.close();

  std::vector<std::unique_ptr<ExperimentStrategy>> strategies;
  strategies.push_back(std::make_unique<TraditionalExperiment>());
  strategies.push_back(std::make_unique<SimpleParallelExperiment>());
  strategies.push_back(std::make_unique<AdvancedParallelExperiment>());

  for (int s = 0; s < n_samples; ++s) {
    std::cout << "Starting Sample " << s + 1 << "/" << n_samples << "..."
              << std::endl;
    auto cbn = CBN::cbn_generator(topology, n_networks, n_vars, 1, 1);

    ExperimentResults last_results;

    for (const auto &strat : strategies) {
      last_results = strat->run(cbn);
      log_to_csv(output_file, last_results, s + 1, topology, n_networks,
                 n_vars);

      std::cout << "  [Strategy Processed] " << last_results.strategy_name
                << " in " << last_results.total_ms << " ms" << std::endl;
    }
    persist_sample(s + 1, cbn, last_results);
  }
  return 0;
}