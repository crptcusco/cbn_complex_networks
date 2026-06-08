#include <iostream>
#include "cbnetwork/campaign_manager.hpp"

using namespace cbnetwork;

int main(int argc, char** argv) {
    std::string config_file = "campaign_config.json";
    if (argc > 1) {
        config_file = argv[1];
    }

    std::cout << "Autonomous Campaign Execution System" << std::endl;
    std::cout << "Loading config from: " << config_file << std::endl;

    auto experiments = CampaignManager::load_campaign_config(config_file);
    if (experiments.empty()) {
        std::cerr << "No experiments found in config or file could not be read." << std::endl;
        return 1;
    }

    CampaignManager::run_campaign(experiments);

    std::cout << "Campaign Execution Finished." << std::endl;
    return 0;
}
