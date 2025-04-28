#pragma once
#include <vector>
using MatrixD = std::vector<std::vector<double>>;

// solve one LP tuple
double solve_lp_for_tuple(
    const MatrixD& G,
    int a_idx, int b_idx,
    const std::vector<int>& X,
    const std::vector<int>& psi,
    bool debug=false);

// full m-height via enumeration (parallelized)
double compute_m_height(
    const MatrixD& G,
    int m,
    bool debug=false);
