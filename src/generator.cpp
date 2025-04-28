#include "generator.h"

MatrixD construct_generator_matrix(const MatrixI &P) {
    int k = P.size();
    int n_minus_k = P[0].size();
    int n = k + n_minus_k;
    MatrixD G(k, std::vector<double>(n, 0.0));
    for (int i = 0; i < k; ++i) G[i][i] = 1.0;
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < n_minus_k; ++j)
            G[i][k + j] = static_cast<double>(P[i][j]);
    return G;
}
