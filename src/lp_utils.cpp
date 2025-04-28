// src/lp_utils.cpp

#include "lp_utils.h"
#include "utils.h"
#include <glpk.h>
#include <set>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>
#include <omp.h>

double solve_lp_for_tuple(
  const MatrixD& G,
  int a_idx, int b_idx,
  const std::vector<int>& X,
  const std::vector<int>& psi,
  bool debug)
{
    // silence unused‐param warning
    (void)debug;

    int k = G.size(), n = (int)G[0].size();
    // build Y = all columns \ {a,b} \ X
    std::set<int> all_cols;
    for (int i = 0; i < n; ++i) all_cols.insert(i);
    all_cols.erase(a_idx);
    all_cols.erase(b_idx);
    for (int j : X) all_cols.erase(j);
    std::vector<int> Y(all_cols.begin(), all_cols.end());

    // turn off GLPK terminal output
    glp_term_out(GLP_OFF);

    // create LP
    glp_prob *lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);

    // k free variables u0..u(k-1)
    glp_add_cols(lp, k);
    for (int i = 1; i <= k; ++i) {
        glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, i, psi[0] * G[i-1][a_idx]);
    }

    // total constraints = 2*|X| + 2 (norm) + 2*|Y|
    int total_rows = 2 * (int)X.size() + 2 + 2 * (int)Y.size();
    glp_add_rows(lp, total_rows);

    std::vector<int>    idx; idx.reserve(k+1);
    std::vector<double> val; val.reserve(k+1);
    int row = 1;

    // --- X constraints ---
    for (int j : X) {
        int r = std::find(X.begin(), X.end(), j) - X.begin() + 1;
        //  (ψ[r]*G[*][j] - ψ[0]*G[*][a]) · u  ≤ 0
        glp_set_row_bnds(lp, row, GLP_UP, 0.0, 0.0);
        idx = {0}; val = {0.0};
        for (int i = 0; i < k; ++i) {
            idx.push_back(i+1);
            val.push_back(psi[r]*G[i][j] - psi[0]*G[i][a_idx]);
        }
        glp_set_mat_row(lp, row, (int)idx.size()-1, idx.data(), val.data());
        ++row;

        //  -(ψ[r]*G[*][j]) · u  ≤ -1
        glp_set_row_bnds(lp, row, GLP_UP, 0.0, -1.0);
        idx = {0}; val = {0.0};
        for (int i = 0; i < k; ++i) {
            idx.push_back(i+1);
            val.push_back(-psi[r]*G[i][j]);
        }
        glp_set_mat_row(lp, row, (int)idx.size()-1, idx.data(), val.data());
        ++row;
    }

    // --- normalization on b_idx: sum G[*][b]·u = 1 via two ineqs ---
    //  sum G[*][b]·u ≤ 1
    glp_set_row_bnds(lp, row, GLP_UP, 0.0, 1.0);
    idx = {0}; val = {0.0};
    for (int i = 0; i < k; ++i) {
        idx.push_back(i+1);
        val.push_back(G[i][b_idx]);
    }
    glp_set_mat_row(lp, row, (int)idx.size()-1, idx.data(), val.data());
    ++row;
    //  sum G[*][b]·u ≥ 1
    glp_set_row_bnds(lp, row, GLP_LO, 1.0, 0.0);
    idx = {0}; val = {0.0};
    for (int i = 0; i < k; ++i) {
        idx.push_back(i+1);
        val.push_back(G[i][b_idx]);
    }
    glp_set_mat_row(lp, row, (int)idx.size()-1, idx.data(), val.data());
    ++row;

    // --- Y constraints: ± sum G[*][j]·u ≤ 1 ---
    for (int j : Y) {
        // +sum ≤ 1
        glp_set_row_bnds(lp, row, GLP_UP, 0.0, 1.0);
        idx = {0}; val = {0.0};
        for (int i = 0; i < k; ++i) {
            idx.push_back(i+1);
            val.push_back(G[i][j]);
        }
        glp_set_mat_row(lp, row, (int)idx.size()-1, idx.data(), val.data());
        ++row;
        // -sum ≤ 1
        glp_set_row_bnds(lp, row, GLP_UP, 0.0, 1.0);
        idx = {0}; val = {0.0};
        for (int i = 0; i < k; ++i) {
            idx.push_back(i+1);
            val.push_back(-G[i][j]);
        }
        glp_set_mat_row(lp, row, (int)idx.size()-1, idx.data(), val.data());
        ++row;
    }

    // solve
    glp_simplex(lp, nullptr);
    int status = glp_get_status(lp);
    double z;
    if      (status == GLP_NOFEAS) z = 0.0;
    else if (status == GLP_UNBND ) z = std::numeric_limits<double>::infinity();
    else                            z = glp_get_obj_val(lp);

    glp_delete_prob(lp);
    return z;
}

double compute_m_height(
  const MatrixD& G,
  int m,
  bool debug)
{
    int n = (int)G[0].size();
    double max_v = -std::numeric_limits<double>::infinity();
    #pragma omp parallel for collapse(2) reduction(max:max_v) schedule(dynamic)
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        if (a == b) continue;
        std::vector<int> cols(n);
        std::iota(cols.begin(), cols.end(), 0);
        cols.erase(std::remove(cols.begin(), cols.end(), a), cols.end());
        cols.erase(std::remove(cols.begin(), cols.end(), b), cols.end());

        std::vector<bool> mask(cols.size());
        std::fill(mask.begin(), mask.begin() + (m-1), true);
        do {
          std::vector<int> X;
          for (size_t i = 0; i < cols.size(); ++i)
            if (mask[i]) X.push_back(cols[i]);

          for (int pm = 0; pm < (1 << (m-1)); ++pm) {
            std::vector<int> psi(m,1);
            for (int bit = 0; bit < m-1; ++bit)
              psi[bit+1] = (pm & (1<<bit)) ? 1 : -1;
            double z = solve_lp_for_tuple(G, a, b, X, psi, false);
            if (std::isinf(z)) {
              max_v = z;  // inf will dominate the reduction
            } else {
              max_v = std::max(max_v, z);
            }
          }
        } while (std::prev_permutation(mask.begin(), mask.end()));
      }
    }

    if (debug) std::cerr << "m-height = " << max_v << "\n";
    return max_v < 0 ? 0.0 : max_v;
}
