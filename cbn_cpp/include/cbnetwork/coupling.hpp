#ifndef COUPLING_HPP
#define COUPLING_HPP

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

namespace cbnetwork {

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

} // namespace cbnetwork

#endif // COUPLING_HPP
