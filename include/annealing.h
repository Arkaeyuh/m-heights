#pragma once
#include <vector>
#include <random>
#include <utility>
using MatrixI = std::vector<std::vector<int>>;

// standard annealing + circulant‐initial helpers
std::pair<MatrixI,double> simulated_annealing_run(
    const MatrixI &P0, int m,
    double T0, double alpha,
    int iter_per_temp, double Tmin,
    bool debug=false);

MatrixI generate_initial_P(int k, int nk, std::mt19937 &rng);
MatrixI generate_circulant_P(int k, int nk, std::mt19937 &rng);
MatrixI generate_broken_circulant_P(int k, int nk, std::mt19937 &rng);
MatrixI generate_noisy_circulant_P(int k, int nk, std::mt19937 &rng);
MatrixI generate_broken_circulant_P2(int k, int nk, std::mt19937 &rng);