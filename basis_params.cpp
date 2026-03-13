#include "basis_params.hpp"
#include <random>
#include <cmath>

namespace ecg1d {

BasisParams BasisParams::conj_params() const {
    BasisParams result;
    result.u = std::conj(u);
    result.A = A.conjugate();
    result.B = B.conjugate();
    result.R = R.conjugate();
    result.name = name;
    return result;
}

BasisParams BasisParams::from_arrays(Cd u, const MatrixXcd& A, const MatrixXcd& B,
                                     const VectorXcd& R, int name) {
    BasisParams bp;
    bp.u = u;
    bp.A = A;
    bp.B = B;
    bp.R = R;
    bp.name = name;
    return bp;
}

// Helper: generate standard normal samples matching numpy's default_rng(seed).standard_normal()
// NumPy uses PCG64 + ziggurat. We cannot reproduce it exactly in C++.
// For testing, we'll load reference values from Python. For production, use mt19937.
static void fill_normal(std::mt19937_64& gen, double* data, int n) {
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n; i++) {
        data[i] = dist(gen);
    }
}

BasisParams BasisParams::randoman(int N, std::optional<int> seed, int name,
                                  double margin, double take_diag_fraction) {
    unsigned long long s = seed.has_value() ? static_cast<unsigned long long>(*seed) : 42;
    std::mt19937_64 gen(s);
    std::normal_distribution<double> dist(0.0, 1.0);

    // Random complex (N,N) matrix X
    MatrixXcd X(N, N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            X(i, j) = Cd(dist(gen), dist(gen));

    // Random diagonal d
    VectorXcd d(N);
    for (int i = 0; i < N; i++)
        d(i) = Cd(dist(gen), dist(gen));

    MatrixXcd BBB = MatrixXcd::Zero(N, N);
    for (int i = 0; i < N; i++)
        BBB(i, i) = d(i);

    // H0 = (X + X^T)/2 + BBB  (symmetric, not Hermitian)
    MatrixXcd H0 = (X + X.transpose()) / 2.0 + BBB;

    // Compute minimum eigenvalue of Hermitian part
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H0);
    double w_min = es.eigenvalues()(0);
    double alpha = (margin - w_min) + 0.1;
    if (alpha < 0) alpha = 0.1;

    MatrixXcd K = H0 + alpha * MatrixXcd::Identity(N, N);

    // Re-check
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es2(K);
    double w_min_K = es2.eigenvalues()(0);
    if (w_min_K < margin) {
        K = K + (margin - w_min_K + 1e-12) * MatrixXcd::Identity(N, N);
    }

    // B = diag(take_diag_fraction * diag(K))
    MatrixXcd B = MatrixXcd::Zero(N, N);
    for (int i = 0; i < N; i++)
        B(i, i) = take_diag_fraction * K(i, i);

    MatrixXcd A = K - B;

    // R
    VectorXcd R(N);
    for (int i = 0; i < N; i++)
        R(i) = Cd(dist(gen), dist(gen));

    // u
    Cd u_val(dist(gen), dist(gen));

    return from_arrays(u_val, A, B, R, name);
}

} // namespace ecg1d
