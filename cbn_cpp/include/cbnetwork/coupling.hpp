#ifndef COUPLING_HPP
#define COUPLING_HPP

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

namespace cbnetwork {

inline void get_combinations(int n, int k, std::vector<std::vector<int>>& out) {
    std::vector<int> p(k);
    std::iota(p.begin(), p.end(), 0);
    while (true) {
        out.push_back(p);
        int i = k - 1;
        while (i >= 0 && p[i] == n - k + i) i--;
        if (i < 0) break;
        p[i]++;
        for (int j = i + 1; j < k; j++) p[j] = p[i] + j - i;
    }
}

class CouplingStrategy {
public:
    virtual ~CouplingStrategy() = default;
    virtual std::string generate_coupling_function(const std::vector<int>& output_variables) const = 0;
    virtual std::vector<std::vector<int>> to_cnf(const std::vector<int>& output_variables, int coupling_variable) const = 0;
};

class OrCoupling : public CouplingStrategy {
public:
    std::string generate_coupling_function(const std::vector<int>& output_variables) const override {
        if (output_variables.empty()) return " ";
        std::string result = " ";
        for (size_t i = 0; i < output_variables.size(); ++i) {
            result += std::to_string(output_variables[i]);
            if (i < output_variables.size() - 1) result += " ∨ ";
        }
        result += " ";
        return result;
    }

    std::vector<std::vector<int>> to_cnf(const std::vector<int>& output_variables, int coupling_variable) const override {
        std::vector<std::vector<int>> clauses;
        for (int var : output_variables) {
            clauses.push_back({-var, coupling_variable});
        }
        std::vector<int> last_clause = output_variables;
        last_clause.push_back(-coupling_variable);
        clauses.push_back(last_clause);
        return clauses;
    }
};

class AndCoupling : public CouplingStrategy {
public:
    std::string generate_coupling_function(const std::vector<int>& output_variables) const override {
        if (output_variables.empty()) return " ";
        std::string result = " ";
        for (size_t i = 0; i < output_variables.size(); ++i) {
            result += std::to_string(output_variables[i]);
            if (i < output_variables.size() - 1) result += " ∧ ";
        }
        result += " ";
        return result;
    }

    std::vector<std::vector<int>> to_cnf(const std::vector<int>& output_variables, int coupling_variable) const override {
        std::vector<std::vector<int>> clauses;
        for (int var : output_variables) {
            clauses.push_back({var, -coupling_variable});
        }
        std::vector<int> last_clause;
        for (int var : output_variables) {
            last_clause.push_back(-var);
        }
        last_clause.push_back(coupling_variable);
        clauses.push_back(last_clause);
        return clauses;
    }
};

class ThresholdCoupling : public CouplingStrategy {
public:
    int threshold;
    ThresholdCoupling(int k) : threshold(k) {}

    std::string generate_coupling_function(const std::vector<int>& output_variables) const override {
        std::string result = " Threshold(" + std::to_string(threshold) + ", {";
        for (size_t i = 0; i < output_variables.size(); ++i) {
            result += std::to_string(output_variables[i]);
            if (i < output_variables.size() - 1) result += ", ";
        }
        result += "}) ";
        return result;
    }

    std::vector<std::vector<int>> to_cnf(const std::vector<int>& output_variables, int coupling_variable) const override {
        int n = output_variables.size();
        int k = threshold;
        std::vector<std::vector<int>> clauses;

        if (k > n) {
            clauses.push_back({-coupling_variable});
            return clauses;
        }

        // Implication 1: (sum(Vi) >= k) => C
        std::vector<std::vector<int>> combinations_k;
        get_combinations(n, k, combinations_k);
        for (const auto& combo : combinations_k) {
            std::vector<int> clause;
            for (int idx : combo) {
                clause.push_back(-output_variables[idx]);
            }
            clause.push_back(coupling_variable);
            clauses.push_back(clause);
        }

        // Implication 2: C => (sum(Vi) >= k)
        std::vector<std::vector<int>> combinations_nk;
        get_combinations(n, n - k + 1, combinations_nk);
        for (const auto& combo : combinations_nk) {
            std::vector<int> clause;
            clause.push_back(-coupling_variable);
            for (int idx : combo) {
                clause.push_back(output_variables[idx]);
            }
            clauses.push_back(clause);
        }

        return clauses;
    }
};

} // namespace cbnetwork

#endif // COUPLING_HPP
