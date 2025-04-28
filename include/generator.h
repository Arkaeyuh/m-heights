#ifndef GENERATOR_H
#define GENERATOR_H

#include <vector>
using MatrixI = std::vector<std::vector<int>>;
using MatrixD = std::vector<std::vector<double>>;

// Build G = [I_k | P]
MatrixD construct_generator_matrix(const MatrixI &P);

#endif // GENERATOR_H
