# Kick Convergence Test — Session Log (2026-04-08)

## Goal

Implement a convergence test that sweeps `T_pulse` values to verify that TDVP-based finite-duration pulse kicks converge to the exact instantaneous kick result as `T_pulse → 0`.

## Changes Made

### 1. `kicking_scale` in `HamiltonianTerms` (src/tdvp_solver.hpp)

Added `double kicking_scale = 1.0` to `HamiltonianTerms`. This multiplies the kicking potential, allowing `V = (kappa/T_pulse) * cos(2*k_L*z)` via `kicking_scale = 1/T_pulse`.

Applied in three places:
- `compute_total_energy()` in `src/tdvp_solver.cpp`
- `grad_H_for_alpha()` in `src/tdvp_solver.cpp`
- `compute_HS_ij()` in `src/svm.cpp`

Default value `1.0` preserves all existing behavior.

### 2. Bug fix: u-coefficient evolution in real-time TDVP (src/tdvp_solver.cpp)

**Bug**: `update_basis_function()` had `if (idx.a1 == 1 || ...) continue;` which unconditionally skipped u parameters. This meant real-time TDVP (used in `realtime_tdvp_step`) never updated the linear combination coefficients u.

**Fix**: Removed the `idx.a1 == 1` check. For imaginary-time TDVP, u entries are already skipped via `dz_idx < 0` (because `updata_constant = K`). For real-time TDVP (`updata_constant = 0`), u is now correctly updated.

### 3. Convergence test function (main.cpp)

`run_kick_convergence_test(int n_kicks)` implements:

| Step | Description |
|------|-------------|
| A | Build N=1 ground state: K=5 Gaussians with A ∈ {0.15, 0.5, 1.5, 5.0, 15.0}, + imaginary-time TDVP polish |
| B | Exact reference via `kicked_exact_1particle()` |
| C | Direct multiplication baseline: `apply_analytic_kick` + `free_evolve_fixed_basis` |
| D | Sweep T_pulse ∈ {0.5, 0.2, 0.1, 0.05, 0.02, 0.01}: TDVP evolution through pulse |
| E | Print convergence table with relative errors and convergence order |

### 4. CLI flag (main.cpp)

`--kick-convergence` triggers `run_kick_convergence_test(n_kicks)`. Uses `--n-kicks` to set kick count (default 50).

## Files Modified

| File | Change |
|------|--------|
| `src/tdvp_solver.hpp` | Added `kicking_scale` field to `HamiltonianTerms` |
| `src/tdvp_solver.cpp` | Applied `kicking_scale` in 2 functions; fixed u update in `update_basis_function` |
| `src/svm.cpp` | Applied `kicking_scale` in `compute_HS_ij` |
| `main.cpp` | Added `run_kick_convergence_test()` + `--kick-convergence` CLI flag |

## Current Results

The test compiles and runs. Ground state is exact (E = 0.5, error ≈ 3e-15). Direct multiplication baseline gives reasonable results (rel_error ≈ 12% after 3 kicks).

**However, the TDVP pulse results are numerically unstable** — energies diverge instead of converging. This is likely due to:

1. **Stiff dynamics**: As `T_pulse → 0`, `kicking_scale = 1/T_pulse` → ∞, creating an extremely strong potential that overwhelms the RK2 midpoint integrator
2. **Basis limitation**: K=5 pure Gaussians (B=0, R=0 initially) may lack variational freedom for momentum-carrying states
3. **Step size**: `dt_kick = T_pulse/20` might be too coarse for the stiff regime

## Next Steps

Possible approaches to stabilize the TDVP pulse convergence:

1. **Higher-order integrator**: Replace RK2 midpoint with RK4 or implicit integrator for stiff dynamics
2. **Adaptive step control**: Reduce dt automatically when energy drift exceeds threshold
3. **Richer basis**: Start with momentum-augmented basis (nonzero B, R) so the TDVP has more variational freedom from the start
4. **Moderate T_pulse regime**: Focus on T_pulse ≥ 0.1 where the potential is not as extreme, and increase basis size
5. **Symplectic integrator**: The TDVP equations have symplectic structure; exploit it for better long-time stability
