#include "generator.h"
#include "lp_utils.h"
#include "annealing.h"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

int main() {
    int n = 10, k = 5, m = 5;
    int nk = n - k;
    std::mt19937 rng(std::random_device{}());

    int tries = 0;
    MatrixI P0;
    double h0;

    // keep drawing random P0 until its m-height is finite
    do {
    P0 = generate_circulate_P(k, nk, rng);
    MatrixD G0 = construct_generator_matrix(P0);
    h0 = compute_m_height(G0, m, /* debug = */ true);
    ++tries;
    std::cout << "m-height=" << h0 << "\n";
    std::cout << vec_to_string(P0) << "\n";
    } while (std::isinf(h0) || h0 > 350);

    // report the first finite one
    std::cout << "Found finite initial m-height=" << h0 
            << " after " << tries << " tries\n";

    // now run annealing from that P0
    auto [P_best, h_best] = simulated_annealing_run(
        P0, m, /*T0=*/50.0, /*alpha=*/0.9, /*iters=*/10, /*Tmin=*/1.0, /*debug=*/true
    );

    std::cout << "Best m-height = " << h_best << "\n";
    std::cout << "P matrix:\n";
    for (auto &row : P_best) {
        for (int v : row) std::cout << v << ' ';
        std::cout << "\n";
    }

    std::ofstream ofs("results.json");
    ofs << "{\n  \"mHeight\": " << h_best << ",\n  \"P\": [\n";
    for (int i = 0; i < k; ++i) {
        ofs << "    [";
        for (int j = 0; j < nk; ++j)
            ofs << P_best[i][j] << (j+1<nk?", ":"");
        ofs << "]" << (i+1<k?",":"") << "\n";
    }
    ofs << "  ]\n}";


    // MatrixI P = {
    //     {  4, -4, -3, -3,  1 },
    //     {  4,  2,  2,  4, -1 },
    //     {  1,  1,  3, -1,  2 },
    //     { -3, -1, -4,  4, -5 },
    //     { -3,  2, -2,  5,  2 }
    // };
    // MatrixD G = construct_generator_matrix(P);
    // std::cout << "m-height = " << compute_m_height(G, 5, /* debug = */ true) << "\n";

}
