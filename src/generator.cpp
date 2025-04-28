#include "generator.h"

MatrixD construct_generator_matrix(const MatrixI &P) {
    int k = P.size();
    int nk = P[0].size();
    int n = k + nk;
    MatrixD G(k, std::vector<double>(n, 0.0));
    // identity part
    for (int i = 0; i < k; ++i) G[i][i] = 1.0;
    // P part
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < nk; ++j)
            G[i][k+j] = double(P[i][j]);
    return G;
}
