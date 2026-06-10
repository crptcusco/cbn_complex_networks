#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/directededge.hpp"
#include "cbnetwork/internalvariable.hpp"
#include "cbnetwork/localnetwork.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

using namespace cbnetwork;

int main() {
  std::freopen("verificacion_debug.log", "w", stdout);

  // 1. Inicialización de Redes
  auto net1 = std::make_shared<LocalNetwork>(1, std::vector<int>{1, 2, 3});
  auto net2 = std::make_shared<LocalNetwork>(2, std::vector<int>{4, 5, 6, 7});
  auto net3 = std::make_shared<LocalNetwork>(3, std::vector<int>{8, 9, 10});
  auto net4 =
      std::make_shared<LocalNetwork>(4, std::vector<int>{11, 12, 13, 14});

  // 2. Definición de Acoplamientos (Relations: (2,1), (3,2), (2,3), (4,3))
  // Nota: El constructor de DirectedEdge suele ser: (input, var_idx, output,
  // local_var_idx, input_vars, logic)
  auto e21 = std::make_shared<DirectedEdge>(2, 15, 1, 1, std::vector<int>{4, 5},
                                            " 4 OR 5 ");
  auto e32 = std::make_shared<DirectedEdge>(3, 16, 2, 4, std::vector<int>{8, 9},
                                            " 8 OR 9 ");
  auto e23 = std::make_shared<DirectedEdge>(2, 17, 3, 8, std::vector<int>{6, 7},
                                            " 6 OR 7 ");
  auto e43 = std::make_shared<DirectedEdge>(
      4, 18, 3, 10, std::vector<int>{13, 14}, " 13 OR 14 ");

  // 3. Registrar señales
  net2->output_signals.push_back(e21);
  net1->process_input_signals({e21});
  net3->output_signals.push_back(e32);
  net2->process_input_signals({e32});
  net2->output_signals.push_back(e23);
  net3->process_input_signals({e23});
  net4->output_signals.push_back(e43);
  net3->process_input_signals({e43});

  std::vector<std::shared_ptr<LocalNetwork>> networks = {net1, net2, net3,
                                                         net4};
  std::vector<std::shared_ptr<DirectedEdge>> edges = {e21, e32, e23, e43};

  // 4. Ejecución del Solver
  CBN cbn(networks, edges);
  cbn.find_local_attractors_sequential();
  cbn.find_compatible_pairs();

  // ... (después de cbn.find_compatible_pairs();)

  std::cout << "\n--- DIAGNOSTICO DE COMPATIBILIDAD ---" << std::endl;
  for (const auto &edge : cbn.l_directed_edges) {
    std::cout << "Relacion entre Net " << edge->input_local_network << " y Net "
              << edge->output_local_network << std::endl;
    // Accedemos a la estructura donde se guardan los pares compatibles
    // Si esto es 0, el problema está en la lógica de las funciones booleanas
    std::cout << "Pares compatibles encontrados: "
              << edge->d_comp_pairs_attractors_by_value.size() << std::endl;
  }

  // --- DIAGNOSTICO DE REDES ---
  for (const auto &net : cbn.l_local_networks) {
    std::cout << "Red " << net->index << " tiene " << net->attractor_count
              << " atractores locales." << std::endl;
  }

  cbn.mount_stable_attractor_fields();

  std::cout << "SISTEMA COMPLETO CARGADO." << std::endl;
  std::cout << "Total de Campos de Atractores encontrados: "
            << cbn.d_attractor_fields.size() << std::endl;

  return 0;
}