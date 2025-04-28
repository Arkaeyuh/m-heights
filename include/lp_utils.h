#ifndef LP_UTILS_H
#define LP_UTILS_H

#include "utils.h"
#include "generator.h"
#include <vector>

// Solve one LP tuple (a,b,X,psi).  Returns objective (∞ if unbounded).
double solve_lp_for_tuple(const MatrixD& G,
                          int a_idx, int b_idx,
                          const std::vector<int>& X,
                          const std::vector<int>& psi,
                          bool debug=false);

// Brute‐force over all tuples to compute the m-height.
double compute_m_height(const MatrixD &G, int m, bool debug=false);

#endif // LP_UTILS_H
