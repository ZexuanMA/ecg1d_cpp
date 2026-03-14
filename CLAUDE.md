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

## Coding Conventions

- **Types**: `Cd` = `complex<double>`, `MatrixXcd`/`VectorXcd` = Eigen types (defined in `types.hpp`)
- **Naming**: snake_case for functions, CamelCase for structs/classes
- **Headers**: `#pragma once`
- **Style**: const-correct, pass by const reference, inline math comments
- **Physical constants**: `hbar=1, mass=1, omega=1` (in `physical_constants.hpp`)
