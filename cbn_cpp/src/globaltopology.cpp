#include "cbnetwork/globaltopology.hpp"
#include <memory>

namespace cbnetwork {

PathDigraph::PathDigraph(int n_nodes) : GlobalTopology(4, {}) {
    for (int i = 1; i < n_nodes; ++i) {
        l_edges.push_back({i, i + 1});
    }
}

CycleDigraph::CycleDigraph(int n_nodes) : GlobalTopology(3, {}) {
    for (int i = 1; i < n_nodes; ++i) {
        l_edges.push_back({i, i + 1});
    }
    l_edges.push_back({n_nodes, 1});
}

CompleteDigraph::CompleteDigraph(int n_nodes) : GlobalTopology(1, {}) {
    for (int i = 1; i <= n_nodes; ++i) {
        for (int j = 1; j <= n_nodes; ++j) {
            if (i != j) {
                l_edges.push_back({i, j});
            }
        }
    }
}

std::shared_ptr<GlobalTopology> GlobalTopology::generate_sample_topology(int v_topology, int n_nodes) {
    if (v_topology == 1) return std::make_shared<CompleteDigraph>(n_nodes);
    if (v_topology == 3) return std::make_shared<CycleDigraph>(n_nodes);
    if (v_topology == 4) return std::make_shared<PathDigraph>(n_nodes);
    return nullptr;
}

} // namespace cbnetwork
