#include "generator.h"
#include "lp_utils.h"
#include "annealing.h"
#include "utils.h"
#include <iostream>
#include <random>
#include <fstream>
#include <cmath>
#include <omp.h>

void local_optimize(MatrixI &P, int m){
    bool improved=true;
    while(improved){
      improved=false;
      auto G = construct_generator_matrix(P);
      double base = compute_m_height(G,m,false);
      // try flipping each entry by ±1
      for(size_t i=0;i<P.size();++i) for(size_t j=0;j<P[0].size();++j){
        for(int d : {+1,-1}){
          P[i][j] += d;
          double h = compute_m_height(construct_generator_matrix(P),m,false);
          if(h < base){
            base=h; improved=true;
          } else {
            P[i][j] -= d;
          }
        }
        if(improved) break;
      }
    }
}

int main(){
    const int n = 9, k = 5, m = 4;
    const int nk = n - k;
    std::mt19937 rng(std::random_device{}());

    // 1) Find a good circulant starting P0
    MatrixI P0;
    double   h0;
    int tries = 0;
    // do {
    //     P0 = generate_broken_circulant_P2(k, nk, rng);
    //     h0 = compute_m_height(
    //         construct_generator_matrix(P0),
    //         m,
    //         /*debug=*/true
    //     );
    //     ++tries;
    // } while (std::isinf(h0) || h0 > 100);  // only accept reasonably small

    P0 = {
        { -5,  4, -4,  3 },
        { -3,  5, -1, -4 },
        {  3,  2, -5, -8 },
        {  4, -1, -3, -4 },
        { -3,  3,  1, -3 }
    };
    h0 = compute_m_height(
        construct_generator_matrix(P0),
        m,
        /*debug=*/true
    );

    std::cout << "Found initial finite m-height=" << h0
              << " after " << tries << " tries\n";
    std::cout << "P0 matrix:\n" << vec_to_string(P0) << "\n";

    // 2) Run _num_runs_ parallel annealing chains and pick the best
    int num_runs = omp_get_max_threads();
    MatrixI bestP = P0;
    double   bestH = h0;

#pragma omp parallel
    {
        auto result = simulated_annealing_run(
            P0,
            m,
            /*T0=*/50,
            /*alpha=*/0.95,
            /*iter_per_temp=*/30,
            /*Tmin=*/1,
            /*debug=*/true
        );

    #pragma omp critical
        {
            if (result.second < bestH) {
                bestH = result.second;
                bestP = result.first;
            }
        }
    }

    // 3) Print out final best
    std::cout << "\n==== Final Best after SA ====\n";
    std::cout << "Best m-height = " << bestH << "\n";
    std::cout << "P matrix:\n";
    for (auto &row : bestP) {
        for (int v : row) std::cout << v << " ";
        std::cout << "\n";
    }

    // 3) Run local optimization
    local_optimize(bestP, m);
    bestH = compute_m_height(
        construct_generator_matrix(bestP),
        m,
        /*debug=*/true
    );
    std::cout << "\n==== Final Best after local opt ====\n";
    std::cout << "Best m-height = " << bestH << "\n";
    std::cout << "P matrix:\n";
    for (auto &row : bestP) {
        for (int v : row) std::cout << v << " ";
        std::cout << "\n";
    }

    // 4) Save to JSON
    std::ofstream ofs("results.json");
    ofs << "{\n  \"mHeight\": " << bestH << ",\n  \"P\": [\n";
    for (int i = 0; i < k; ++i) {
        ofs << "    [";
        for (int j = 0; j < nk; ++j)
            ofs << bestP[i][j] << (j+1< nk ? ", " : "");
        ofs << "]" << (i+1<k ? "," : "") << "\n";
    }
    ofs << "  ]\n}\n";

    return 0;
}
