#ifndef ANNEALING_H
#define ANNEALING_H

#include "generator.h"
#include <vector>
#include <random>

// tweak one cell ±step, avoid zero‐column
MatrixI perturb_P(const MatrixI &P, int step_size, std::mt19937 &rng);

// simulated annealing driver
std::pair<MatrixI,double> simulated_annealing_run(
    const MatrixI &P_init,
    int m,
    double T0 = 45.0,
    double alpha = 0.9,
    int iter_per_temp = 20,
    double min_temp = 1.0,
    bool debug = false);

// random start P
MatrixI generate_initial_P(int k, int n_minus_k, int m, std::mt19937 &rng);

// random subcirculate P
MatrixI generate_circulate_P(int k, int nk, std::mt19937 &rng);

#endif // ANNEALING_H
