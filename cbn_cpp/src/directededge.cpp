#include "cbnetwork/directededge.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace cbnetwork {

DirectedEdge::DirectedEdge(int idx, int idx_var, int in_net, int out_net,
                           const std::vector<int>& out_vars, const std::string& coup_func)
    : index(idx), index_variable(idx_var), input_local_network(in_net),
      output_local_network(out_net), l_output_variables(out_vars),
      coupling_function(coup_func), kind_signal(2) {

    d_kind_signal[1] = "RESTRICTED";
    d_kind_signal[2] = "NOT COMPUTE";
    d_kind_signal[3] = "STABLE";
    d_kind_signal[4] = "NOT STABLE";

    d_out_value_to_attractor[0] = {};
    d_out_value_to_attractor[1] = {};
    d_comp_pairs_attractors_by_value[0] = {};
    d_comp_pairs_attractors_by_value[1] = {};

    true_table = process_true_table();
}

void DirectedEdge::show() const {
    std::cout << "Index Edge: " << index << " - Relation: " << output_local_network
              << " -> " << input_local_network << " - Variable: " << index_variable << std::endl;
    std::cout << "Variables: [";
    for (size_t i = 0; i < l_output_variables.size(); ++i) {
        std::cout << l_output_variables[i] << (i == l_output_variables.size() - 1 ? "" : ", ");
    }
    std::cout << "], Coupling Function: " << coupling_function << std::endl;
    std::cout << "Kind signal: " << kind_signal << " - " << d_kind_signal.at(kind_signal) << std::endl;
}

void DirectedEdge::show_short() const {
    std::cout << output_local_network << " , " << input_local_network << std::endl;
}

std::pair<int, int> DirectedEdge::get_edge() const {
    return {output_local_network, input_local_network};
}

std::map<std::string, std::string> DirectedEdge::process_true_table() {
    std::map<std::string, std::string> r_true_table;
    int n = l_output_variables.size();
    if (n == 0) return r_true_table;

    // Simple implementation for basic functions like " v1 ∨ v2 " or " v1 ∧ v2 "
    // Since we are translating, let's try to be a bit more robust if possible.
    // For now, let's assume it's OR or AND as per CouplingStrategy.

    for (int i = 0; i < (1 << n); ++i) {
        std::string key = "";
        std::map<char, bool> env;
        for (int j = 0; j < n; ++j) {
            bool val = (i >> (n - 1 - j)) & 1;
            key += val ? "1" : "0";
            env['A' + j] = val;
        }

        bool result = false;
        if (coupling_function.find("∨") != std::string::npos) {
            result = false;
            for (int j = 0; j < n; ++j) result = result || env['A' + j];
        } else if (coupling_function.find("∧") != std::string::npos) {
            result = true;
            for (int j = 0; j < n; ++j) result = result && env['A' + j];
        } else {
            // Default to OR for now if unknown
            result = false;
            for (int j = 0; j < n; ++j) result = result || env['A' + j];
        }

        r_true_table[key] = result ? "1" : "0";
    }

    return r_true_table;
}

} // namespace cbnetwork
