#pragma once
#include "types.hpp"
#include <vector>

namespace ecg1d {

struct PermutationSet {
    int N;
    int SN;  // N!
    std::vector<VectorXi> sigmas;      // each is length-N permutation
    std::vector<int> signs;            // +1 or -1
    std::vector<MatrixXi> matrices;    // N x N permutation matrices

    static PermutationSet generate(int N);

    // sigma_inverse: find j such that sigmas[p][j] == i
    int sigma_inv(int p, int i) const;
};

} // namespace ecg1d
