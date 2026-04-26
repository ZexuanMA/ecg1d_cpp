#pragma once
#include "step4_realtime_tdvp.hpp"
#include "step4_reference_dynamics.hpp"
#include <Eigen/Dense>
#include <string>

namespace ecg1d {
namespace verify {

// Step-4 cross-check: pairwise dynamics diagnostics over snapshot times.
// Writes step4_crosscheck_N{N}.csv (long format: pair,t,quantity,value).
// Compares ECG vs grid, ECG vs DVR, grid vs DVR. For each pair at each
// snapshot time records:
//   delta_E         = E_a - E_b
//   abs_delta_x     = |<x>_a - <x>_b|
//   abs_delta_p     = |<p>_a - <p>_b|
//   L2_nx           = sqrt(integral (n_a(x) - n_b(x))^2 dx)
//   L2_nk           = sqrt(integral (n_a(k) - n_b(k))^2 dk)
//   bhattacharyya   = integral sqrt(n_a(x) n_b(x)) dx   (in [0, 1])
//   fidelity        = |<psi_a|psi_b>|^2 / (norm_a * norm_b)   (N=1 only)
//
// Densities and momentum distributions are first resampled (linear
// interpolation) onto the caller-supplied user x_grid / k_grid so that L^2
// comparisons are meaningful.
void step4_cross_check(int N,
                       const Step4EcgResult& ecg,
                       const Step4RefResult& grid_res,
                       const Step4RefResult& dvr_res,
                       const Eigen::VectorXd& x_grid,
                       const Eigen::VectorXd& k_grid);

} // namespace verify
} // namespace ecg1d
