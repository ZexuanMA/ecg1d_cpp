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

} // namespace verify
} // namespace ecg1d
