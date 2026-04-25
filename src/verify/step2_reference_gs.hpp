#pragma once
#include "types.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace ecg1d {
namespace verify {

struct RefSolverResult {
    double E;                     // N-particle energy (product state: N * E_sp)
    double E_sp;                  // single-particle ground-state energy
    Eigen::VectorXd x;            // real-space grid
    Eigen::VectorXcd psi;         // normalized single-particle ground-state wavefunction on x
    Eigen::VectorXd n_x;          // n(x) = |psi(x)|^2 (single-particle probability density)
    Eigen::VectorXd k;            // momentum grid (matches caller-provided k_grid)
    Eigen::VectorXd n_k;          // n(k) = |psi_k(k)|^2
};

// Step 2a: finite-difference 3-point grid diagonalization of single-particle H_full.
// Reuses the method at main_simple.cpp:345-409.
RefSolverResult reference_gs_grid(int N, int n_grid, double x_max,
                                  const Eigen::VectorXd& k_grid);

// Step 2b: sinc-DVR with closed-form kinetic matrix.
RefSolverResult reference_gs_dvr(int N, int n_dvr, double x_max,
                                 const Eigen::VectorXd& k_grid);

// Writes the reference solver result (energy + n(x) + n(k)) to CSVs under the given
// method tag ("grid" or "dvr") and particle number N.
void write_reference_csvs(const std::string& method_tag, int N,
                          const RefSolverResult& res);

namespace detail {

// Build the FD 3-point kinetic matrix T_jj = 1/(m dx^2), T_{j,j±1} = -1/(2 m dx^2)
// on the regular grid x.
Eigen::MatrixXd build_T_grid_fd(const Eigen::VectorXd& x);

// Build the sinc-DVR kinetic matrix on the regular grid x.
//   T_ii = pi^2 / (6 m dx^2);  T_ij = (-1)^(i-j) / (m dx^2 (i-j)^2)
Eigen::MatrixXd build_T_dvr_sinc(const Eigen::VectorXd& x);

// V_full(x_j) = 1/2 m omega^2 x_j^2 + kappa * cos(2 k_L x_j)  (Step-1/2 Hamiltonian)
Eigen::VectorXd build_V_full(const Eigen::VectorXd& x);

// V_ho(x_j) = 1/2 m omega^2 x_j^2  (Step-4 evolution: cos quench off)
Eigen::VectorXd build_V_ho_only(const Eigen::VectorXd& x);

// Normalize psi so that sum_j |psi_j|^2 * dx = 1.
void normalize_grid_state(Eigen::VectorXcd& psi, double dx);

// Direct DFT of psi on x onto k_grid; returns n(k) = |psi_tilde(k)|^2.
Eigen::VectorXd compute_nk(const Eigen::VectorXcd& psi,
                           const Eigen::VectorXd& x,
                           const Eigen::VectorXd& k_grid);

} // namespace detail

} // namespace verify
} // namespace ecg1d
