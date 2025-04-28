#include "annealing.h"
#include "lp_utils.h"
#include "utils.h"
#include "generator.h"
#include <cmath>
#include <random>
#include <iostream>

MatrixI perturb_P(const MatrixI &P, int step_size, std::mt19937 &rng) {
    std::uniform_int_distribution<int> dr(0, P.size()-1);
    std::uniform_int_distribution<int> dc(0, P[0].size()-1);
    std::uniform_int_distribution<int> ds(0,1);
    MatrixI Q = P;
    int i = dr(rng), j = dc(rng);
    int δ = ds(rng) ? +step_size : -step_size;
    Q[i][j] += δ;
    // keep columns non‐zero
    bool all_zero = true;
    for (int r = 0; r < (int)Q.size(); ++r)
      if (Q[r][j] != 0) { all_zero = false; break; }
    if (all_zero) Q[i][j] -= δ;
    return Q;
}

std::pair<MatrixI,double> simulated_annealing_run(
    const MatrixI &P0, int m,
    double T0, double alpha,
    int iter_per_temp, double Tmin,
    bool debug)
{
    std::mt19937 rng(std::random_device{}());
    MatrixI cur = P0, best = P0;
    double fcur = compute_m_height(
        construct_generator_matrix(cur), m, debug);
    double fbest = fcur, T = T0;

    if (debug) std::cerr << "init height=" << fcur << "\n";
    if (debug) std::cerr << "init P=" << vec_to_string(cur) << "\n";
    if (debug) std::cerr << "T0=" << T0 << "\n";
    if (debug) std::cerr << "Tmin=" << Tmin << "\n";
    if (debug) std::cerr << "alpha=" << alpha << "\n";
    if (debug) std::cerr << "iter_per_temp=" << iter_per_temp << "\n";

    while (T > Tmin) {
        int accepted = 0;
        for (int it = 0; it < iter_per_temp; ++it) {
            auto cand = perturb_P(cur, 1, rng);
            double fcand = compute_m_height(
                construct_generator_matrix(cand), m, debug);
            double Δ = fcand - fcur;
            if (Δ < 0 || std::uniform_real_distribution<>(0,1)(rng) < std::exp(-Δ/T)) {
                cur = cand; fcur = fcand; ++accepted;
                if (fcand < fbest) {
                    best = cand; fbest = fcand;
                    if (debug) {
                        std::cerr << "new best=" << fbest << "\n";
                        std::cerr << "new P=" << vec_to_string(best) << "\n";
                    }
                }
            }
        }
        if (debug) std::cerr << "T=" << T << " accepted=" << accepted << "\n";
        if (debug) std::cerr << "fcur=" << fcur << "\n";
        if (debug) std::cerr << "cur P=" << vec_to_string(cur) << "\n";
        double rate = double(accepted)/iter_per_temp;
        if      (rate < 0.2) T *= alpha*alpha;
        else if (rate > 0.8) T *= std::sqrt(alpha);
        else                 T *= alpha;
    }
    return {best,fbest};
}

MatrixI generate_initial_P(int k, int nk, std::mt19937 &rng) {
    std::uniform_int_distribution<int> d(-10,10);
    MatrixI P(k, std::vector<int>(nk));
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < nk; ++j)
        P[i][j] = d(rng);
    // no zero‐column
    for (int j = 0; j < nk; ++j) {
      bool all_zero = true;
      for (int i = 0; i < k; ++i)
        if (P[i][j] != 0) { all_zero = false; break; }
      if (all_zero) P[0][j] = 1;
    }
    return P;
}

MatrixI generate_circulant_P(int k, int nk, std::mt19937 &rng) {
    std::uniform_int_distribution<int> d(-30,30);
    std::vector<int> base(nk);
    bool all_zero = true;
    for (int j = 0; j < nk; ++j) {
        base[j] = d(rng);
        if (base[j] != 0) all_zero = false;
    }
    if (all_zero) base[0] = 1;
    MatrixI P(k, std::vector<int>(nk));
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < nk; ++j)
        P[i][j] = base[(j + i) % nk];
    return P;
}

/**
 * Generate a “broken” circulant P of size k×nk:
 *   1) start with a true circulant of width nk,
 *   2) compute g = gcd(k, nk), which is the period,
 *   3) for each residue‐class r=0..g-1, perturb one entry in each
 *      of the subsequent rows of that class by ±1,
 *   4) ensure no column of P is ever identically zero.
 */
MatrixI generate_broken_circulant_P(int k, int nk, std::mt19937 &rng) {
  std::uniform_int_distribution<int> val_dist(-20, 20);
  std::uniform_int_distribution<int> sign_dist(0, 1);
  std::uniform_int_distribution<int> col_dist(0, nk - 1);

  // 1) pick a random base row
  std::vector<int> base(nk);
  bool all_zero = true;
  for (int j = 0; j < nk; ++j) {
      base[j] = val_dist(rng);
      if (base[j] != 0) all_zero = false;
  }
  if (all_zero) base[0] = 1;

  // 2) build true circulant P
  MatrixI P(k, std::vector<int>(nk));
  for (int i = 0; i < k; ++i)
      for (int j = 0; j < nk; ++j)
          P[i][j] = base[(j + i) % nk];

  // 3) break each repeat‐class mod g = gcd(k, nk)
  int g = std::gcd(k, nk);
  for (int r = 0; r < g; ++r) {
      // collect all rows i ≡ r (mod g)
      std::vector<int> rows;
      for (int i = r; i < k; i += g)
          rows.push_back(i);

      // for each *later* row in that class, perturb one entry
      for (size_t idx = 1; idx < rows.size(); ++idx) {
          int irow = rows[idx];
          int jcol = col_dist(rng);
          int delta = sign_dist(rng) ? +1 : -1;
          P[irow][jcol] += delta;
      }
  }

  // 4) ensure no all-zero column
  for (int j = 0; j < nk; ++j) {
      bool col_all_zero = true;
      for (int i = 0; i < k; ++i) {
          if (P[i][j] != 0) {
              col_all_zero = false;
              break;
          }
      }
      if (col_all_zero) {
          // just fix the first row in that column
          P[0][j] = 1;
      }
  }

  return P;
}

MatrixI generate_broken_circulant_P2(int k,int nk, std::mt19937 &rng){
  // like a circulant base, but now each column‐shift you also
  // randomly swap two entries (with low probability) to break symmetry.
  std::uniform_int_distribution<int>   base_dist(-20,20);
  std::bernoulli_distribution          glitch(0.2), swapero(0.1);
  std::vector<int> base(nk);
  bool allz = true;
  for(int j=0;j<nk;++j){
    base[j] = base_dist(rng);
    if(base[j]!=0) allz=false;
  }
  if(allz) base[0]=1;
  MatrixI P(k,std::vector<int>(nk));
  for(int i=0;i<k;++i){
    for(int j=0;j<nk;++j){
      P[i][j] = base[(j+i)%nk];
      if(glitch(rng)) P[i][j] += (rng()%2?1:-1);
    }
    if(swapero(rng)){
      int c1 = rng()%nk, c2 = rng()%nk;
      std::swap(P[i][c1],P[i][c2]);
    }
  }
  // ensure no zero‐col
  for(int j=0;j<nk;++j){
    bool z=true;
    for(int i=0;i<k;++i) if(P[i][j]!=0){ z=false; break; }
    if(z) P[0][j]=1;
  }
  return P;
}

// “Noisy circulant”: base row in [–3..3], then each cell is base[(j+i)%nk]
// plus a ±1 “glitch” with probability p.
MatrixI generate_noisy_circulant_P(int k, int nk, std::mt19937 &rng) {
  std::uniform_int_distribution<int>  base_dist(-20, 20);
  std::bernoulli_distribution         glitch(0.25);     // 25% of entries get perturbed
  std::uniform_int_distribution<int>  sign_dist(0, 1);

  // pick a random base row of length nk
  std::vector<int> base(nk);
  bool all_zero = true;
  for (int j = 0; j < nk; ++j) {
      base[j] = base_dist(rng);
      if (base[j] != 0) all_zero = false;
  }
  if (all_zero) base[0] = 1;  // avoid the trivial all-zero

  // build the “circulant+noise” P
  MatrixI P(k, std::vector<int>(nk));
  for (int i = 0; i < k; ++i) {
      for (int j = 0; j < nk; ++j) {
          int v = base[(j + i) % nk];
          if (glitch(rng)) {
              // randomly bump it by ±1
              v += (sign_dist(rng) ? +1 : -1);
          }
          P[i][j] = v;
      }
  }

  // make sure no column is identically zero
  for (int j = 0; j < nk; ++j) {
      bool col_zero = true;
      for (int i = 0; i < k; ++i) {
          if (P[i][j] != 0) { col_zero = false; break; }
      }
      if (col_zero) P[0][j] = 1;
  }

  return P;
}