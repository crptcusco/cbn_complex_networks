#include "cbnetwork/globaltopology.hpp"
#include <memory>
#include <random>
#include <algorithm>
#include <numeric>
#include <set>

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

AleatoryFixedDigraph::AleatoryFixedDigraph(int n_nodes, int n_edges) : GlobalTopology(2, {}) {
    if (n_nodes <= 1) return;
    std::random_device rd;
    std::mt19937 gen(rd());

    // Spanning tree to ensure weak connectivity
    std::vector<int> nodes(n_nodes);
    std::iota(nodes.begin(), nodes.end(), 1);
    std::shuffle(nodes.begin(), nodes.end(), gen);

    for (int i = 1; i < n_nodes; ++i) {
        std::uniform_int_distribution<> dis(0, i - 1);
        l_edges.push_back({nodes[dis(gen)], nodes[i]});
    }

    // Add remaining edges
    std::set<std::pair<int, int>> edge_set;
    for (auto& e : l_edges) edge_set.insert(e);

    int target_edges = std::min(n_edges, n_nodes * (n_nodes - 1));
    std::uniform_int_distribution<> dis_node(1, n_nodes);

    while ((int)l_edges.size() < target_edges) {
        int u = dis_node(gen);
        int v = dis_node(gen);
        if (u != v && edge_set.find({u, v}) == edge_set.end()) {
            l_edges.push_back({u, v});
            edge_set.insert({u, v});
        }
    }
}

DorogovtsevMendesDigraph::DorogovtsevMendesDigraph(int n_nodes) : GlobalTopology(7, {}) {
    if (n_nodes < 3) return;
    l_edges = {{1, 2}, {2, 3}, {3, 1}};
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 4; i <= n_nodes; ++i) {
        std::uniform_int_distribution<> dis(0, (int)l_edges.size() - 1);
        auto edge = l_edges[dis(gen)];
        l_edges.push_back({i, edge.first});
        l_edges.push_back({i, edge.second});
    }
}

SmallWorldGraph::SmallWorldGraph(int n_nodes, int k_neighbors, double p_rewire) : GlobalTopology(8, {}) {
    if (n_nodes <= k_neighbors) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_prob(0.0, 1.0);
    std::uniform_int_distribution<> dis_node(1, n_nodes);

    std::set<std::pair<int, int>> edge_set;
    for (int i = 1; i <= n_nodes; ++i) {
        for (int j = 1; j <= k_neighbors / 2; ++j) {
            int neighbor = i + j;
            if (neighbor > n_nodes) neighbor -= n_nodes;
            edge_set.insert({std::min(i, neighbor), std::max(i, neighbor)});
        }
    }

    // Rewire
    std::vector<std::pair<int, int>> edges(edge_set.begin(), edge_set.end());
    for (auto& edge : edges) {
        if (dis_prob(gen) < p_rewire) {
            edge_set.erase(edge);
            int u = edge.first;
            int v = dis_node(gen);
            while (v == u || edge_set.count({std::min(u, v), std::max(u, v)})) {
                v = dis_node(gen);
            }
            edge_set.insert({std::min(u, v), std::max(u, v)});
        }
    }

    for (auto& e : edge_set) l_edges.push_back(e);
}

ScaleFreeGraph::ScaleFreeGraph(int n_nodes, int m_edges) : GlobalTopology(9, {}) {
    if (n_nodes <= m_edges) return;
    std::random_device rd;
    std::mt19937 gen(rd());

    // Initial complete graph of m_edges + 1 nodes
    for (int i = 1; i <= m_edges + 1; ++i) {
        for (int j = i + 1; j <= m_edges + 1; ++j) {
            l_edges.push_back({i, j});
        }
    }

    for (int i = m_edges + 2; i <= n_nodes; ++i) {
        std::vector<int> degrees;
        for (int node = 1; node < i; ++node) {
            int deg = 0;
            for (auto& e : l_edges) if (e.first == node || e.second == node) deg++;
            degrees.insert(degrees.end(), deg, node);
        }

        std::set<int> targets;
        while ((int)targets.size() < m_edges) {
            std::uniform_int_distribution<> dis(0, (int)degrees.size() - 1);
            targets.insert(degrees[dis(gen)]);
        }
        for (int t : targets) l_edges.push_back({i, t});
    }
}

RandomGraph::RandomGraph(int n_nodes, double p_edge) : GlobalTopology(10, {}) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 1; i <= n_nodes; ++i) {
        for (int j = i + 1; j <= n_nodes; ++j) {
            if (dis(gen) < p_edge) {
                l_edges.push_back({i, j});
            }
        }
    }
}

std::shared_ptr<GlobalTopology> GlobalTopology::generate_sample_topology(int v_topology, int n_nodes, int n_edges) {
    if (v_topology == 1) return std::make_shared<CompleteDigraph>(n_nodes);
    if (v_topology == 2) return std::make_shared<AleatoryFixedDigraph>(n_nodes, n_edges == -1 ? n_nodes : n_edges);
    if (v_topology == 3) return std::make_shared<CycleDigraph>(n_nodes);
    if (v_topology == 4) return std::make_shared<PathDigraph>(n_nodes);
    if (v_topology == 7) return std::make_shared<DorogovtsevMendesDigraph>(n_nodes);
    if (v_topology == 8) return std::make_shared<SmallWorldGraph>(n_nodes, 3, 0.5);
    if (v_topology == 9) return std::make_shared<ScaleFreeGraph>(n_nodes, 2);
    if (v_topology == 10) return std::make_shared<RandomGraph>(n_nodes, 0.5);
    return nullptr;
}

} // namespace cbnetwork
