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

    static std::shared_ptr<GlobalTopology> generate_sample_topology(int v_topology, int n_nodes);
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

} // namespace cbnetwork

#endif // GLOBALTOPOLOGY_HPP
