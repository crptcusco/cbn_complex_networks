#ifndef GLOBALTOPOLOGY_HPP
#define GLOBALTOPOLOGY_HPP

#include <vector>
#include <utility>
#include <memory>

namespace cbnetwork {

class GlobalTopology {
public:
    int v_topology;
    std::vector<std::pair<int, int>> l_edges;

    GlobalTopology(int v_topo, const std::vector<std::pair<int, int>>& edges)
        : v_topology(v_topo), l_edges(edges) {}

    virtual ~GlobalTopology() = default;

    static std::shared_ptr<GlobalTopology> generate_sample_topology(int v_topology, int n_nodes, int n_edges = -1);
};

class PathDigraph : public GlobalTopology {
public:
    PathDigraph(int n_nodes);
};

class CycleDigraph : public GlobalTopology {
public:
    CycleDigraph(int n_nodes);
};

class CompleteDigraph : public GlobalTopology {
public:
    CompleteDigraph(int n_nodes);
};

class AleatoryFixedDigraph : public GlobalTopology {
public:
    AleatoryFixedDigraph(int n_nodes, int n_edges);
};

class DorogovtsevMendesDigraph : public GlobalTopology {
public:
    DorogovtsevMendesDigraph(int n_nodes);
};

class SmallWorldGraph : public GlobalTopology {
public:
    SmallWorldGraph(int n_nodes, int k_neighbors, double p_rewire);
};

class ScaleFreeGraph : public GlobalTopology {
public:
    ScaleFreeGraph(int n_nodes, int m_edges);
};

class RandomGraph : public GlobalTopology {
public:
    RandomGraph(int n_nodes, double p_edge);
};

} // namespace cbnetwork

#endif // GLOBALTOPOLOGY_HPP
