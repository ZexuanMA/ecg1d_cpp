#pragma once
#include "step1_imag_time_ecg.hpp"
#include "step2_reference_gs.hpp"
#include <Eigen/Dense>

namespace ecg1d {
namespace verify {

// Step 3: cross-check ECG vs reference ground states.
// Writes step3_crosscheck_N{N}.csv and prints a summary table.
// Requires Step 1 and Step 2 already in memory; no disk re-reads.
//
// ecg_psi_x: optional complex ECG wavefunction on ecg_x_grid (size 0 -> skip the
// complex-overlap column; only meaningful for N=1).
void step3_cross_check(int N,
                       const Step1Result& ecg_res,
                       const Eigen::VectorXd& ecg_x_grid,
                       const Eigen::VectorXd& ecg_n_x,
                       const Eigen::VectorXd& ecg_k_grid,
                       const Eigen::VectorXd& ecg_n_k,
                       const Eigen::VectorXcd& ecg_psi_x,
                       const RefSolverResult& grid_res,
                       const RefSolverResult& dvr_res);

} // namespace verify
} // namespace ecg1d
