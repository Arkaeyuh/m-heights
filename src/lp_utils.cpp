#include "lp_utils.h"
#include "utils.h"
#include <glpk.h>
#include <set>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>

double solve_lp_for_tuple(const MatrixD& G,
  int a_idx, int b_idx,
  const std::vector<int>& X,
  const std::vector<int>& psi,
  bool debug)
{
  int k = G.size(), n = G[0].size();
  // build Y = all columns minus {a,b} minus X
  std::set<int> all_cols;
  for (int i = 0; i < n; ++i) all_cols.insert(i);
  all_cols.erase(a_idx);
  all_cols.erase(b_idx);
  for (int j : X) all_cols.erase(j);
  std::vector<int> Y(all_cols.begin(), all_cols.end());

  // turn off GLPK terminal output
  glp_term_out(GLP_OFF);

  glp_prob* lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  // k free variables u0..u(k-1)
  glp_add_cols(lp, k);
  for (int i = 1; i <= k; ++i) {
      glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
      glp_set_obj_coef(lp, i, psi[0] * G[i-1][a_idx]);
  }

  // total constraints = 2*|X| + 2 (normalization) + 2*|Y|
  int total_rows = 2*X.size() + 2 + 2*Y.size();
  glp_add_rows(lp, total_rows);

  std::vector<int>    idx; idx.reserve(k+1);
  std::vector<double> val; val.reserve(k+1);
  int row = 1;

  // --- X constraints ---
  for (int j : X) {
      int r = std::find(X.begin(), X.end(), j) - X.begin() + 1;

      // (ψ[r]*G[*][j] - ψ[0]*G[*][a]) · u   ≤ 0
      glp_set_row_bnds(lp, row, GLP_UP, /*lb*/0.0, /*ub*/0.0);
      idx = {0}; val = {0.0};
      for (int i = 0; i < k; ++i) {
          idx.push_back(i+1);
          val.push_back(psi[r]*G[i][j] - psi[0]*G[i][a_idx]);
      }
      glp_set_mat_row(lp, row, (int)idx.size() - 1, idx.data(), val.data());
      ++row;

      // - (ψ[r]*G[*][j]) · u  ≤ -1
      glp_set_row_bnds(lp, row, GLP_UP, /*lb*/0.0, /*ub*/-1.0);
      idx = {0}; val = {0.0};
      for (int i = 0; i < k; ++i) {
          idx.push_back(i+1);
          val.push_back(-psi[r]*G[i][j]);
      }
      glp_set_mat_row(lp, row, (int)idx.size() - 1, idx.data(), val.data());
      ++row;
  }

  // --- normalization on b_idx:  sum_i G[i][b] u_i = 1  via two inequalities ---

  //  sum_i G[*][b] · u   ≤ 1
  glp_set_row_bnds(lp, row, GLP_UP, /*lb*/0.0, /*ub*/1.0);
  idx = {0}; val = {0.0};
  for (int i = 0; i < k; ++i) {
      idx.push_back(i+1);
      val.push_back(G[i][b_idx]);
  }
  glp_set_mat_row(lp, row, (int)idx.size() - 1, idx.data(), val.data());
  ++row;

  //  sum_i G[*][b] · u   ≥ 1
  glp_set_row_bnds(lp, row, GLP_LO, /*lb*/1.0, /*ub*/0.0);
  idx = {0}; val = {0.0};
  for (int i = 0; i < k; ++i) {
      idx.push_back(i+1);
      val.push_back(G[i][b_idx]);
  }
  // note: for GLP_LO the sign stays positive, so no need to negate
  glp_set_mat_row(lp, row, (int)idx.size() - 1, idx.data(), val.data());
  ++row;

  // --- Y constraints:  ± sum_i G[i][j] u_i  ≤ 1  ---
  for (int j : Y) {
      //  sum_i G[*][j] · u   ≤ 1
      glp_set_row_bnds(lp, row, GLP_UP, /*lb*/0.0, /*ub*/1.0);
      idx = {0}; val = {0.0};
      for (int i = 0; i < k; ++i) {
          idx.push_back(i+1);
          val.push_back(G[i][j]);
      }
      glp_set_mat_row(lp, row, (int)idx.size() - 1, idx.data(), val.data());
      ++row;

      // - sum_i G[*][j] · u   ≤ 1
      glp_set_row_bnds(lp, row, GLP_UP, /*lb*/0.0, /*ub*/1.0);
      idx = {0}; val = {0.0};
      for (int i = 0; i < k; ++i) {
          idx.push_back(i+1);
          val.push_back(-G[i][j]);
      }
      glp_set_mat_row(lp, row, (int)idx.size() - 1, idx.data(), val.data());
      ++row;
  }

  // solve
  glp_simplex(lp, nullptr);

  // check status properly (==, not =!)
  int st = glp_get_status(lp);
  double z;
  if      (st == GLP_NOFEAS) z = 0.0;
  else if (st == GLP_UNBND ) z = std::numeric_limits<double>::infinity();
  else                        z = glp_get_obj_val(lp);

  glp_delete_prob(lp);
  return z;
}


double compute_m_height(const MatrixD &G, int m, bool debug) {
    int k = G.size(), n = G[0].size();
    double max_v = -std::numeric_limits<double>::infinity();
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        if (a==b) continue;
        // build cols = [0..n)
        std::vector<int> cols(n);
        std::iota(cols.begin(), cols.end(), 0);
        cols.erase(std::remove(cols.begin(), cols.end(), a), cols.end());
        cols.erase(std::remove(cols.begin(), cols.end(), b), cols.end());

        // choose m-1 of them
        std::vector<bool> mask(cols.size());
        std::fill(mask.begin(), mask.begin()+m-1, true);
        do {
          std::vector<int> X;
          for (size_t i=0; i<cols.size(); ++i)
            if (mask[i]) X.push_back(cols[i]);

          for (int pm=0; pm< (1<<(m-1)); ++pm) {
            std::vector<int> psi(m,1);
            for (int bit=0; bit<m-1; ++bit)
              psi[bit+1] = (pm&(1<<bit))?1:-1;
            double z = solve_lp_for_tuple(G,a,b,X,psi,debug);
            if (std::isinf(z)) return std::numeric_limits<double>::infinity();
            max_v = std::max(max_v, z);
          }
        } while(std::prev_permutation(mask.begin(), mask.end()));
      }
    }
    if(debug) {
        std::cerr << "m-height: " << max_v << "\n";
    }
    return max_v < 0 ? 0.0 : max_v;
}
