#ifndef NETWORK_FACTORY_HPP
#define NETWORK_FACTORY_HPP

#include <memory>
#include "cbnetwork/cbnetwork.hpp"
#include "cbnetwork/campaign_manager.hpp"

namespace cbnetwork {

class NetworkFactory {
public:
    static std::shared_ptr<CBN> create_cbn(const ExperimentConfig& config);
};

} // namespace cbnetwork

#endif // NETWORK_FACTORY_HPP
