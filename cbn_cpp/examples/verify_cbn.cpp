#include <iostream>
#include <vector>
#include <memory>
#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/localnetwork.hpp"
#include "cbnetwork/directededge.hpp"
#include "cbnetwork/internalvariable.hpp"

using namespace cbnetwork;

int main() {
    // A very simple 2-network system
    // Network 1: x1(t+1) = x1(t)
    // Network 2: x2(t+1) = x1(t) (coupled)

    // Network 1
    auto net1 = std::make_shared<LocalNetwork>(1, std::vector<int>{1});
    auto var1 = std::make_shared<InternalVariable>(1, std::vector<std::vector<int>>{{1}});
    net1->descriptive_function_variables.push_back(var1);

    // Network 2
    auto net2 = std::make_shared<LocalNetwork>(2, std::vector<int>{2});
    // For net2, x2(t+1) = x3(t), where x3 is the signal from net1
    auto var2 = std::make_shared<InternalVariable>(2, std::vector<std::vector<int>>{{3}});
    net2->descriptive_function_variables.push_back(var2);

    // Directed Edge from Net 1 to Net 2
    // Output variable from Net 1 is x1
    // Signal variable index is 3
    auto edge12 = std::make_shared<DirectedEdge>(1, 3, 2, 1, std::vector<int>{1}, " 1 ");

    net1->process_input_signals({});
    net1->output_signals.push_back(edge12);
    net2->process_input_signals({edge12});

    std::vector<std::shared_ptr<LocalNetwork>> networks = {net1, net2};
    std::vector<std::shared_ptr<DirectedEdge>> edges = {edge12};

    CBN cbn(networks, edges);

    std::cout << "Finding local attractors..." << std::endl;
    cbn.find_local_attractors_sequential();
    cbn.show_local_attractors();

    std::cout << "Finding compatible pairs..." << std::endl;
    cbn.find_compatible_pairs();
    cbn.show_attractor_pairs();

    std::cout << "Mounting attractor fields..." << std::endl;
    cbn.mount_stable_attractor_fields();
    cbn.show_stable_attractor_fields();

    return 0;
}
