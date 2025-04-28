#include "annealing.h"
#include "lp_utils.h"
#include "utils.h"
#include <cmath>
#include <random>
#include <iostream>

MatrixI perturb_P(const MatrixI &P, int step_size, std::mt19937 &rng) {
    std::uniform_int_distribution<int> dr(0, P.size()-1);
    std::uniform_int_distribution<int> dc(0, P[0].size()-1);
    std::uniform_int_distribution<int> ds(0,1);
    MatrixI Q = P;
    int i = dr(rng), j = dc(rng);
    int δ = ds(rng)? +step_size : -step_size;
    Q[i][j] += δ;
    bool zero = true;
    for (auto &row : Q) if (row[j]!=0) { zero=false; break; }
    if (zero) Q[i][j] -= δ; 
    return Q;
}

std::pair<MatrixI,double> simulated_annealing_run(
    const MatrixI &P0, int m,
    double T0, double α, int iters, double Tmin, bool dbg)
{
    std::mt19937 rng(std::random_device{}());
    MatrixI cur = P0, best = P0;
    double fcur = compute_m_height(construct_generator_matrix(cur),m,dbg);
    double fbest = fcur;
    double T = T0;
    if (dbg) std::cerr<<"init height="<<fcur<<"\n";
    if(dbg) std::cerr<<"init P="<<vec_to_string(cur)<<"\n";
    while (T > Tmin) {
      int acc=0;
      for (int i=0; i<iters; ++i) {
        auto cand = perturb_P(cur,1,rng);
        double fcand = compute_m_height(construct_generator_matrix(cand),m,dbg);
        double Δ = fcand - fcur;
        std::uniform_real_distribution<> u(0,1);
        if (Δ<0 || u(rng)<std::exp(-Δ/T)) {
          cur = cand; fcur = fcand; ++acc;
          if (fcur < fbest) { best=cur; fbest=fcur;
            if (dbg) std::cerr<<"new best="<<fbest<<"\n";
            if (dbg) std::cerr<<"new P="<<vec_to_string(best)<<"\n";
          }
        }
      }
      double r = double(acc)/iters;
      if (r<0.2)     T *= α*α;
      else if (r>0.8)T *= std::sqrt(α);
      else           T *= α;
    }
    return {best,fbest};
}

MatrixI generate_initial_P(int k, int nk, int m, std::mt19937 &rng) {
    std::uniform_int_distribution<int> d(-3,3);
    MatrixI P(k, std::vector<int>(nk));
    for(int i=0;i<k;++i)for(int j=0;j<nk;++j)P[i][j]=d(rng);
    for(int j=0;j<nk;++j){
      bool zero=true;
      for(int i=0;i<k;++i) if(P[i][j]!=0){zero=false;break;}
      if(zero)P[0][j]=1;
    }
    return P;
}

MatrixI generate_circulate_P(int k, int nk, std::mt19937 &rng) {
  // pick a random base row in [-3 ,3]
  std::uniform_int_distribution<int> dist(-3,3);
  std::vector <int> base(nk);
  bool all_zero = true;
  for (int j = 0; j < nk; ++j) {
    base[j] = dist(rng);
    if (base[j] != 0) all_zero = false;
  }
  // avoid all-zero case
  if(all_zero) {
    base[0] = 1;
  }

  // build circulant matrix

  MatrixI P(k, std::vector<int>(nk));
  for(int i = 0; i < k; ++i) {
    for (int j = 0; j < nk; ++j) {
      P[i][j] = base[(j + i) % nk];
    }
  }

  return P;
  
}