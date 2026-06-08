#ifndef CAMPAIGN_MANAGER_HPP
#define CAMPAIGN_MANAGER_HPP

#include <string>
#include <vector>
#include <memory>
#include "cbnetwork/cbnetwork.hpp"

namespace cbnetwork {

struct ExperimentConfig {
    std::string experiment_id;
    int topology_type;
    int num_networks;
    int vars_per_network;
    int input_variables;
    int output_variables;
    int max_clauses;
    int max_literals;
    std::string coupling_strategy_name;
    int threshold_value = 0;
};

class CampaignManager {
public:
    static std::vector<ExperimentConfig> load_campaign_config(const std::string& filepath);
    static void run_campaign(const std::vector<ExperimentConfig>& experiments);
    static void export_experiment_data(const ExperimentConfig& config,
                                     std::shared_ptr<CBN> cbn,
                                     double execution_time_s);
};

} // namespace cbnetwork

#endif // CAMPAIGN_MANAGER_HPP
