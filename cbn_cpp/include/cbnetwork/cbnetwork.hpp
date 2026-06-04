#ifndef CBN_HPP
#define CBN_HPP

#include <vector>
#include <memory>
#include <map>
#include <set>
#include <string>
#include <tuple>
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
    std::map<int, std::shared_ptr<LocalAttractor>> d_local_attractors_ptr;
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
    bool update_network_by_index(std::shared_ptr<LocalNetwork> o_local_network_update);

    void find_local_attractors_sequential();
    void find_local_attractors_parallel();
    void find_local_attractors_parallel_with_weights();
    void find_local_attractors_brute_force_turbo();

    void find_compatible_pairs();
    void find_compatible_pairs_parallel();
    void find_compatible_pairs_parallel_with_weights();
    void find_compatible_pairs_turbo();

    void order_edges_by_compatibility();
    void order_edges_by_grade();
    void disorder_edges();

    void mount_stable_attractor_fields();
    void mount_stable_attractor_fields_parallel();
    void mount_stable_attractor_fields_parallel_chunks();
    void mount_stable_attractor_fields_turbo();

    void generate_attractor_dictionary();
    void process_kind_signal(std::shared_ptr<LocalNetwork> o_local_network);
    void generate_global_scenes();
    void count_fields_by_global_scenes();

    std::vector<std::shared_ptr<LocalAttractor>> get_attractors_by_input_signal_value(int index_variable_signal, int signal_value);
    std::shared_ptr<LocalAttractor> get_local_attractor_by_index(int i_attractor);
    std::shared_ptr<LocalNetwork> get_network_by_index(int index);

    int get_n_local_attractors() const;
    int get_n_pair_attractors() const;
    int get_n_attractor_fields() const;
    int get_n_local_variables() const;

    void show_local_attractors() const;
    void show_attractor_pairs() const;
    void show_stable_attractor_fields() const;
    void show_directed_edges() const;
    void show_coupled_signals_kind() const;
    void show_description() const;
    void show_global_scenes() const;
    void show_resume() const;
    void show_local_attractors_dictionary() const;
    void show_stable_attractor_fields_detailed() const;
    void show_attractor_fields() const;

    void save_attractor_fields_to_json(const std::string& filepath);

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

    // Turbo Kernels
    static std::vector<std::vector<int8_t>> filter_compatible_pairs_kernel(
        const std::vector<std::vector<int>>& fields,
        const std::vector<std::vector<int>>& candidates,
        const std::vector<int>& attr_to_network
    );
};

} // namespace cbnetwork

#endif // CBN_HPP
