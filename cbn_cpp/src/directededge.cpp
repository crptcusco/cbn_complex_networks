#include "cbnetwork/directededge.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stack>
#include <sstream>
#include <memory>

namespace cbnetwork {

class BooleanParser {
public:
    struct Node {
        enum Type { VAR, CONST, NOT, AND, OR, IMPLY, EQUIV } type;
        char var_name;
        bool const_val;
        std::unique_ptr<Node> left, right;
        Node(Type t) : type(t), var_name(0), const_val(false), left(nullptr), right(nullptr) {}
    };

    BooleanParser(const std::string& expr) : s(expr), pos(0) {
        // Replace UTF-8 operators with single characters for easier parsing
        std::string replaced;
        for (size_t i = 0; i < s.length(); ++i) {
            if ((unsigned char)s[i] == 0xE2 && i + 2 < s.length()) {
                if ((unsigned char)s[i+1] == 0x88 && (unsigned char)s[i+2] == 0xA7) { // ∧
                    replaced += '&'; i += 2;
                } else if ((unsigned char)s[i+1] == 0x88 && (unsigned char)s[i+2] == 0xA8) { // ∨
                    replaced += '|'; i += 2;
                } else if ((unsigned char)s[i+1] == 0x86 && (unsigned char)s[i+2] == 0x92) { // →
                    replaced += '>'; i += 2;
                } else if ((unsigned char)s[i+1] == 0x86 && (unsigned char)s[i+2] == 0x94) { // ↔
                    replaced += '='; i += 2;
                } else replaced += s[i];
            } else replaced += s[i];
        }
        s = replaced;
    }

    std::unique_ptr<Node> parse() {
        pos = 0;
        return parse_equivalence();
    }

private:
    std::string s;
    size_t pos;

    void skip_ws() {
        while (pos < s.length() && isspace(s[pos])) pos++;
    }

    // Precedence (lowest to highest): EQUIV (=), IMPLY (>), OR (|), AND (&), NOT (~), Term
    std::unique_ptr<Node> parse_equivalence() {
        auto left = parse_implication();
        skip_ws();
        if (pos < s.length() && s[pos] == '=') {
            pos++;
            auto res = std::make_unique<Node>(Node::EQUIV);
            res->left = std::move(left);
            res->right = parse_equivalence();
            return res;
        }
        return left;
    }

    std::unique_ptr<Node> parse_implication() {
        auto left = parse_disjunction();
        skip_ws();
        if (pos < s.length() && s[pos] == '>') {
            pos++;
            auto res = std::make_unique<Node>(Node::IMPLY);
            res->left = std::move(left);
            res->right = parse_implication();
            return res;
        }
        return left;
    }

    std::unique_ptr<Node> parse_disjunction() {
        auto left = parse_conjunction();
        skip_ws();
        if (pos < s.length() && s[pos] == '|') {
            pos++;
            auto res = std::make_unique<Node>(Node::OR);
            res->left = std::move(left);
            res->right = parse_disjunction();
            return res;
        }
        return left;
    }

    std::unique_ptr<Node> parse_conjunction() {
        auto left = parse_unary();
        skip_ws();
        if (pos < s.length() && s[pos] == '&') {
            pos++;
            auto res = std::make_unique<Node>(Node::AND);
            res->left = std::move(left);
            res->right = parse_conjunction();
            return res;
        }
        return left;
    }

    std::unique_ptr<Node> parse_unary() {
        skip_ws();
        if (pos < s.length() && s[pos] == '~') {
            pos++;
            auto res = std::make_unique<Node>(Node::NOT);
            res->left = parse_unary();
            return res;
        }
        return parse_term();
    }

    std::unique_ptr<Node> parse_term() {
        skip_ws();
        if (pos >= s.length()) return nullptr;
        if (s[pos] == '(') {
            pos++;
            auto res = parse_equivalence();
            skip_ws();
            if (pos < s.length() && s[pos] == ')') pos++;
            return res;
        }
        if (s[pos] == '0' || s[pos] == '1') {
            auto res = std::make_unique<Node>(Node::CONST);
            res->const_val = (s[pos] == '1');
            pos++;
            return res;
        }
        if (isalpha(s[pos])) {
            auto res = std::make_unique<Node>(Node::VAR);
            res->var_name = s[pos];
            pos++;
            return res;
        }
        return nullptr;
    }
};

static bool evaluate_tree(const std::unique_ptr<BooleanParser::Node>& node, const std::map<char, bool>& env) {
    if (!node) return false;
    switch (node->type) {
        case BooleanParser::Node::VAR: return env.at(node->var_name);
        case BooleanParser::Node::CONST: return node->const_val;
        case BooleanParser::Node::NOT: return !evaluate_tree(node->left, env);
        case BooleanParser::Node::AND: return evaluate_tree(node->left, env) && evaluate_tree(node->right, env);
        case BooleanParser::Node::OR: return evaluate_tree(node->left, env) || evaluate_tree(node->right, env);
        case BooleanParser::Node::IMPLY: return !evaluate_tree(node->left, env) || evaluate_tree(node->right, env);
        case BooleanParser::Node::EQUIV: return evaluate_tree(node->left, env) == evaluate_tree(node->right, env);
    }
    return false;
}

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

    std::map<std::string, char> dict_aux_var_saida;
    for (int i = 0; i < n; ++i) {
        dict_aux_var_saida[" " + std::to_string(l_output_variables[i]) + " "] = 'A' + i;
    }

    std::string aux_coupling_function = coupling_function;
    // Replace variable indices with A, B, C...
    // Note: Python replaces " v_idx " with "LETTER"
    for (auto const& [key, val] : dict_aux_var_saida) {
        size_t pos = 0;
        while ((pos = aux_coupling_function.find(key, pos)) != std::string::npos) {
            aux_coupling_function.replace(pos, key.length(), std::string(1, val));
            pos += 1;
        }
    }

    BooleanParser parser(aux_coupling_function);
    auto tree = parser.parse();

    for (int i = 0; i < (1 << n); ++i) {
        std::string key = "";
        std::map<char, bool> env;
        for (int j = 0; j < n; ++j) {
            bool val = (i >> (n - 1 - j)) & 1;
            key += val ? "1" : "0";
            env['A' + j] = val;
        }

        bool result = evaluate_tree(tree, env);
        r_true_table[key] = result ? "1" : "0";
    }

    return r_true_table;
}

} // namespace cbnetwork
