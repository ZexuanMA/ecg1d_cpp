#pragma once
#include "basis_params.hpp"
#include <Eigen/Dense>
#include <vector>

namespace ecg1d {
namespace verify {

// Step 4: real-time TDVP evolution of the Step-1 ground state under H_evolve = T + V_ho
// (cos potential turned off). Writes:
//   - step4_trace_N{N}_K{K}.csv    : t, E, norm, x2, p2, |<psi0|psi_t>|^2
//   - step4_snapshots_N{N}_K{K}.csv: long format, t, kind (x|k), grid, value
// Prints max relative energy drift.
void step4_realtime_tdvp(int N, int K,
                         const std::vector<BasisParams>& basis_init,
                         double T_total, double dt,
                         const Eigen::VectorXd& x_grid,
                         const Eigen::VectorXd& k_grid,
                         int n_snapshots = 5);

} // namespace verify
} // namespace ecg1d
