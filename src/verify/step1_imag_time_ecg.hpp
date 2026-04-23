#pragma once
#include "basis_params.hpp"
#include "tdvp_solver.hpp"
#include <Eigen/Dense>
#include <vector>

namespace ecg1d {
namespace verify {

struct Step1Result {
    std::vector<BasisParams> basis;   // optimized ECG basis
    double E_ecg;                     // final energy
    double E_kinetic;
    double E_harmonic;
    double E_kicking;
    int    n_steps;
};

// Step 1: imaginary-time TDVP ground state of H_full = T + V_ho + kappa*cos(2 k_L x).
// Writes four CSVs (step1_ecg_N{N}_K{K}.csv, step1_ecg_basis_*, step1_ecg_density_*,
// step1_ecg_nk_*) into out/verify/.
Step1Result step1_imag_time_ecg(int N, int K,
                                const Eigen::VectorXd& x_grid,
                                const Eigen::VectorXd& k_grid);

} // namespace verify
} // namespace ecg1d
