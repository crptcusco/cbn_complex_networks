#include "cbnetwork/network_factory.hpp"
#include "cbnetwork/coupling.hpp"

namespace cbnetwork {

std::shared_ptr<CBN> NetworkFactory::create_cbn(const ExperimentConfig& config) {
    std::shared_ptr<CouplingStrategy> strategy;

    if (config.coupling_strategy_name == "OR") {
        strategy = std::make_shared<OrCoupling>();
    } else if (config.coupling_strategy_name == "AND") {
        strategy = std::make_shared<AndCoupling>();
    } else if (config.coupling_strategy_name == "XOR") {
        strategy = std::make_shared<XorCoupling>();
    } else if (config.coupling_strategy_name == "CONST0") {
        strategy = std::make_shared<ConstantCoupling>(0);
    } else if (config.coupling_strategy_name == "CONST1") {
        strategy = std::make_shared<ConstantCoupling>(1);
    } else if (config.coupling_strategy_name == "THRESHOLD") {
        strategy = std::make_shared<ThresholdCoupling>(config.threshold_value);
    } else {
        strategy = std::make_shared<OrCoupling>();
    }

    return CBN::cbn_generator(
        config.topology_type,
        config.num_networks,
        config.vars_per_network,
        config.input_variables,
        config.output_variables,
        config.max_clauses,
        config.max_literals,
        -1,
        strategy
    );
}

} // namespace cbnetwork
