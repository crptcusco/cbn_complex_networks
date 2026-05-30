#ifndef CBN_HPP
#define CBN_HPP

#include <vector>
#include <memory>
#include <map>
#include <set>
#include "cbnetwork/localnetwork.hpp"
#include "cbnetwork/directededge.hpp"
#include "cbnetwork/globaltopology.hpp"
#include "cbnetwork/globalscene.hpp"

namespace cbnetwork {

class CBN {
public:
    std::vector<std::shared_ptr<LocalNetwork>> l_local_networks;
    std::vector<std::shared_ptr<DirectedEdge>> l_directed_edges;
    std::map<int, std::tuple<int, int, int>> d_local_attractors;
    std::map<int, std::vector<int>> d_attractor_fields;
    std::vector<std::shared_ptr<GlobalScene>> l_global_scenes;
    std::map<std::string, int> d_global_scenes_count;
    std::shared_ptr<GlobalTopology> o_global_topology;

    CBN(const std::vector<std::shared_ptr<LocalNetwork>>& networks,
        const std::vector<std::shared_ptr<DirectedEdge>>& edges)
        : l_local_networks(networks), l_directed_edges(edges), o_global_topology(nullptr) {
        process_output_signals();
    }

    void process_output_signals();
    void order_edges_by_grade();
    void disorder_edges();
    void generate_global_scenes();
    void count_fields_by_global_scenes();
    void find_local_attractors_sequential();
    void find_local_attractors_parallel();
    void find_compatible_pairs();
    void find_compatible_pairs_parallel();
    void mount_stable_attractor_fields();
    void mount_stable_attractor_fields_parallel();

    void generate_attractor_dictionary();
    void process_kind_signal(std::shared_ptr<LocalNetwork> o_local_network);

    std::vector<std::shared_ptr<LocalAttractor>> get_attractors_by_input_signal_value(int index_variable_signal, int signal_value);
    std::shared_ptr<LocalAttractor> get_local_attractor_by_index(int i_attractor);

    void show_local_attractors() const;
    void show_attractor_pairs() const;
    void show_stable_attractor_fields() const;

    static std::shared_ptr<CBN> cbn_generator(
        int v_topology,
        int n_local_networks,
        int n_vars_network,
        int n_input_variables,
        int n_output_variables,
        int n_max_of_clauses = 2,
        int n_max_of_literals = 3);

    void _assign_global_indices_to_attractors();
    std::vector<std::string> _generate_local_scenes(std::shared_ptr<LocalNetwork> o_local_network);

private:
    void order_edges_by_compatibility();
};

} // namespace cbnetwork

#endif // CBN_HPP
