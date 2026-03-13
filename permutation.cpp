#include "permutation.hpp"
#include <algorithm>
#include <numeric>

namespace ecg1d {

static int parity(const std::vector<int>& perm) {
    int n = static_cast<int>(perm.size());
    std::vector<bool> visited(n, false);
    int par = 1;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            int cycle_len = 0;
            int j = i;
            while (!visited[j]) {
                visited[j] = true;
                j = perm[j];
                cycle_len++;
            }
            if (cycle_len % 2 == 0) par *= -1;
        }
    }
    return par;
}

PermutationSet PermutationSet::generate(int N) {
    PermutationSet ps;
    ps.N = N;

    // Generate all permutations of [0, 1, ..., N-1]
    std::vector<int> base(N);
    std::iota(base.begin(), base.end(), 0);

    std::vector<std::vector<int>> all_perms;
    do {
        all_perms.push_back(base);
    } while (std::next_permutation(base.begin(), base.end()));

    ps.SN = static_cast<int>(all_perms.size());
    ps.sigmas.resize(ps.SN);
    ps.signs.resize(ps.SN);
    ps.matrices.resize(ps.SN);

    MatrixXi I = MatrixXi::Identity(N, N);

    for (int p = 0; p < ps.SN; p++) {
        // Store sigma
        ps.sigmas[p] = VectorXi(N);
        for (int i = 0; i < N; i++)
            ps.sigmas[p](i) = all_perms[p][i];

        // Sign: currently hardcoded to 1 to match Python's build_sign_dict
        // which returns 1 for all permutations
        ps.signs[p] = 1;

        // Permutation matrix: T[p] = I[sigma[p]]
        // T[p](i, j) = I(sigma[i], j) => row i of T is row sigma[i] of I
        ps.matrices[p] = MatrixXi(N, N);
        for (int i = 0; i < N; i++)
            ps.matrices[p].row(i) = I.row(all_perms[p][i]);
    }

    return ps;
}

int PermutationSet::sigma_inv(int p, int i) const {
    for (int j = 0; j < N; j++) {
        if (sigmas[p](j) == i) return j;
    }
    return -1; // should never happen
}

} // namespace ecg1d
