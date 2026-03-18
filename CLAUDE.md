# CLAUDE.md — ecg1d_cpp

## Project Overview

C++ implementation of ECG1D (Effective Core Gaussian, 1-Dimensional) for quantum many-body simulations of ultracold atom systems. Uses a Gaussian variational ansatz with TDVP (Time-Dependent Variational Principle) evolution.

Companion Python reference implementation lives at `../ecg1d/`.

## Build

```bash
cd build && cmake .. && make
```

- CMake 3.16+, C++17
- Single executable: `build/ecg1d`
- Dependency: Eigen 3.4.0 (auto-fetched via FetchContent)
- Flags: `-O2 -march=native`

## Running

```bash
./build/ecg1d --csv basis_file.csv N K   # Phase 1 & 2 tests (functionals + derivatives)
./build/ecg1d --tdvp                      # 1-particle harmonic TDVP
./build/ecg1d --delta                     # 2-particle delta TDVP
./build/ecg1d --gaussian                  # 2-particle gaussian TDVP
./build/ecg1d --kicking                   # 2-particle kicking TDVP
```

## Validation

Python scripts validate C++ output against the reference implementation:

```bash
python validate.py          # Phase 1: overlap & Hamiltonian functionals
python validate_phase2.py   # Phase 2: variational derivatives
python validate_phase3.py   # Phase 3: gradients & TDVP evolution
```

## Architecture

Three-phase modular design, all in namespace `ecg1d`:

1. **Functionals** (`hamiltonian.hpp/cpp`): Overlap, kinetic, harmonic, delta, gaussian, kicking
2. **Derivatives** (`derivatives.hpp/cpp`, `observable_derivatives.hpp/cpp`): Log-overlap derivatives, `partial_z_first/second`, C matrix
3. **TDVP** (`tdvp_solver.hpp/cpp`, `hamiltonian_gradient.hpp/cpp`): Time evolution via least-squares TDVP step

Key data structures:
- `BasisParams` (`basis_params.hpp`): Variational basis (u, A, B, R, name)
- `PairCache` (`pair_cache.hpp`): Cached pair computations
- `PermutationSet` (`permutation.hpp`): N! permutations with sign factors

## Physics Goal: Quantitative Simulation of Many-Body Dynamical Localization (MBDL)

Reference paper: Guo et al., "Observation of many-body dynamical localization", Science 389, 716 (2025). DOI: 10.1126/science.adn8625

### What the paper does
- Experimentally observes MBDL in 1D Bose gases (Cs atoms) modeled by the Lieb-Liniger QKR (quantum kicked rotor)
- Tunes interaction strength γ from 0 (noninteracting) to 11 (Tonks-Girardeau regime)
- Three experimental signatures of MBDL: (1) momentum distribution n(k) freezing, (2) energy and entropy saturation, (3) distinct decay of first-order correlation function G⁽¹⁾(z)
- Theory modeling captures results only **qualitatively** (quantum Monte Carlo + Floquet methods)

### What we aim to do
- Provide **quantitative** simulation of MBDL using ECG + TDVP — a fundamentally different theoretical approach
- Our code already has the key ingredients: δ contact interaction, kicking potential, TDVP evolution, Gaussian variational basis
- These map directly to the paper's Hamiltonian (Eq. 1): kinetic + periodic kick + external trap + δ interaction

### Strategy
1. Benchmark: N=2 with kicking + δ interaction, validate against exact results
2. Scale up to N=3,4,5,6, observe MBDL signatures (energy saturation, n(k) freezing)
3. Quantitative results for small N would already be a meaningful contribution beyond the paper's qualitative modeling

### Key challenges
- Particle number scaling: N! permutation sum grows fast; N=18 (experiment) is out of reach, but N=3~6 is valuable
- Long-time stability: thousands of TDVP kicks require careful error accumulation control
- Basis completeness: need sufficient K (basis functions) for convergence without computational explosion

## Phase A: Static Optimizer (completed)

See `PHASE_A_RESULTS.md` for full details. Key changes:
- `SolverConfig` with energy convergence, stagnation recovery, adaptive lambda
- `evolution()` dtao_max = 100x initial, diagnostic summary
- `tdvp_step()` always re-solves u, exploratory tiny step on line search failure
- Result: 1-particle error 2.2e-6 → 1.8e-12; Gaussian 8.1e-6 → 4.7e-6; Delta unchanged at 3e-3 (basis limitation)
- Stochastic refine disabled due to numerical instability (accepts ghost energies)

## Coding Conventions

- **Types**: `Cd` = `complex<double>`, `MatrixXcd`/`VectorXcd` = Eigen types (defined in `types.hpp`)
- **Naming**: snake_case for functions, CamelCase for structs/classes
- **Headers**: `#pragma once`
- **Style**: const-correct, pass by const reference, inline math comments
- **Physical constants**: `hbar=1, mass=1, omega=1` (in `physical_constants.hpp`)
