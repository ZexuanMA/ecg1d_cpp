#!/usr/bin/env python3
"""
Generate reference data for C++ validation.
Dumps basis function parameters to CSV and computes all functional values.
"""
import sys
import os
import numpy as np

# Add Python source directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ecg1d'))

from Parameters_z import BasisParams
from Generate_basis_function import generate_basis_function
from overlap import overlap
from Hamiltonian_Kinetic_energy_functional import kinetic_energy_functional
from Hamiltonian_Harmonic_functional import Harmonic_functional
from Hamiltonian_Delta_contact_interation_functional import Delta_contact_interation_functional
from Hamiltonian_Gaussian_interaction_functional import Gaussian_interaction_functional
from Hamiltonian_Kicking_term_functional import kicking_term_functional


def dump_basis_csv(basis_function, filename):
    """Dump basis parameters to CSV: one complex number per line as 'real,imag'"""
    with open(filename, 'w') as f:
        for bp in basis_function:
            f.write(f"{bp.u.real},{bp.u.imag}\n")
            for i in range(bp.N):
                for j in range(bp.N):
                    f.write(f"{bp.A[i,j].real},{bp.A[i,j].imag}\n")
            for i in range(bp.N):
                for j in range(bp.N):
                    f.write(f"{bp.B[i,j].real},{bp.B[i,j].imag}\n")
            for i in range(bp.N):
                f.write(f"{bp.R[i].real},{bp.R[i].imag}\n")
            f.write(f"{bp.name}\n")


def main():
    for N in [1, 2, 3]:
        for basis_n in [1, 2, 3]:
            print(f"\n=== N={N}, K={basis_n} ===")
            basis = generate_basis_function(basis_n, N)

            # Dump parameters
            csv_file = f"basis_N{N}_K{basis_n}.csv"
            dump_basis_csv(basis, csv_file)
            print(f"  Written: {csv_file}")

            # Compute functionals
            S = overlap(basis)
            T = kinetic_energy_functional(basis)
            V = Harmonic_functional(basis)

            print(f"  Overlap:   {S.real:.18e} + {S.imag:.18e}i")
            print(f"  Kinetic:   {T.real:.18e} + {T.imag:.18e}i")
            print(f"  Harmonic:  {V.real:.18e} + {V.imag:.18e}i")

            if N >= 2:
                D = Delta_contact_interation_functional(basis)
                G = Gaussian_interaction_functional(basis)
                print(f"  Delta:     {D.real:.18e} + {D.imag:.18e}i")
                print(f"  Gaussian:  {G.real:.18e} + {G.imag:.18e}i")

            K_val = kicking_term_functional(basis)
            print(f"  Kicking:   {K_val.real:.18e} + {K_val.imag:.18e}i")

            E = (T + V) / S
            print(f"  E(T+V)/S:  {E.real:.18e} + {E.imag:.18e}i")


if __name__ == '__main__':
    main()
