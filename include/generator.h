#pragma once
#include <vector>
using MatrixI = std::vector<std::vector<int>>;
using MatrixD = std::vector<std::vector<double>>;

// build [ I_k | P ] from P
MatrixD construct_generator_matrix(const MatrixI &P);
