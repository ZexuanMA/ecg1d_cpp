#!/usr/bin/env python3
"""
Generate reference data for Phase 2 (derivatives) validation.
Dumps partial_z_first, partial_z_second, and calculate_C values.
"""
import sys
import os
from pathlib import Path
import numpy as np

# Add sibling Python ECG directory to path (robust to custom HOME settings)
ecg_py_dir = Path(__file__).resolve().parent.parent / 'ecg1d'
sys.path.insert(0, str(ecg_py_dir))

from Parameters_z import BasisParams
from Generate_basis_function import generate_basis_function
from Partial_z_first import partial_z_first
from Partial_z_second import partial_z_second
from Calculate_C import calculate_C
from overlap import overlap


def fmt(z):
    return f"{z.real:.18e} + {z.imag:.18e}i"


def main():
    for N in [1, 2]:
        basis_n = 2
        print(f"\n=== N={N}, K={basis_n} ===")
        basis = generate_basis_function(basis_n, N)

        S = overlap(basis)
        print(f"Overlap: {fmt(S)}")

        # Test partial_z_first for various parameter types
        # alpha_1=1 (u derivative), alpha_2=name of basis[0]
        name0 = basis[0].name
        name1 = basis[1].name

        # u derivative (Real=True)
        v = partial_z_first(alpha_1=1, Real=True, basis_function=basis,
                            alpha_2=name0, alpha_3=0, alpha_4=0)
        print(f"pz_first(1,T,{name0},0,0): {fmt(v)}")

        # u derivative (Real=False)
        v = partial_z_first(alpha_1=1, Real=False, basis_function=basis,
                            alpha_2=name0, alpha_3=0, alpha_4=0)
        print(f"pz_first(1,F,{name0},0,0): {fmt(v)}")

        # B derivative (alpha_1=2)
        v = partial_z_first(alpha_1=2, Real=True, basis_function=basis,
                            alpha_2=name0, alpha_3=0, alpha_4=0)
        print(f"pz_first(2,T,{name0},0,0): {fmt(v)}")

        v = partial_z_first(alpha_1=2, Real=False, basis_function=basis,
                            alpha_2=name0, alpha_3=0, alpha_4=0)
        print(f"pz_first(2,F,{name0},0,0): {fmt(v)}")

        # R derivative (alpha_1=3)
        v = partial_z_first(alpha_1=3, Real=True, basis_function=basis,
                            alpha_2=name0, alpha_3=0, alpha_4=0)
        print(f"pz_first(3,T,{name0},0,0): {fmt(v)}")

        v = partial_z_first(alpha_1=3, Real=False, basis_function=basis,
                            alpha_2=name1, alpha_3=0, alpha_4=0)
        print(f"pz_first(3,F,{name1},0,0): {fmt(v)}")

        # A derivative (alpha_1=4, diagonal)
        v = partial_z_first(alpha_1=4, Real=True, basis_function=basis,
                            alpha_2=name0, alpha_3=0, alpha_4=0)
        print(f"pz_first(4,T,{name0},0,0): {fmt(v)}")

        if N >= 2:
            # A derivative (off-diagonal)
            v = partial_z_first(alpha_1=4, Real=True, basis_function=basis,
                                alpha_2=name0, alpha_3=0, alpha_4=1)
            print(f"pz_first(4,T,{name0},0,1): {fmt(v)}")

        # Test partial_z_second
        # B-B second derivative
        v = partial_z_second(alpha_1=2, beta_1=2, basis_function=basis,
                             alpha_2=name0, alpha_3=0, alpha_4=0,
                             beta_2=name0, beta_3=0, beta_4=0)
        print(f"pz_second(2,2,{name0},0,0,{name0},0,0): {fmt(v)}")

        # u-u second derivative
        v = partial_z_second(alpha_1=1, beta_1=1, basis_function=basis,
                             alpha_2=name0, alpha_3=0, alpha_4=0,
                             beta_2=name0, beta_3=0, beta_4=0)
        print(f"pz_second(1,1,{name0},0,0,{name0},0,0): {fmt(v)}")

        # u-B second derivative
        v = partial_z_second(alpha_1=1, beta_1=2, basis_function=basis,
                             alpha_2=name0, alpha_3=0, alpha_4=0,
                             beta_2=name0, beta_3=0, beta_4=0)
        print(f"pz_second(1,2,{name0},0,0,{name0},0,0): {fmt(v)}")

        # R-B second derivative
        v = partial_z_second(alpha_1=3, beta_1=2, basis_function=basis,
                             alpha_2=name0, alpha_3=0, alpha_4=0,
                             beta_2=name0, beta_3=0, beta_4=0)
        print(f"pz_second(3,2,{name0},0,0,{name0},0,0): {fmt(v)}")

        # A-B second derivative (diagonal A)
        v = partial_z_second(alpha_1=4, beta_1=2, basis_function=basis,
                             alpha_2=name0, alpha_3=0, alpha_4=0,
                             beta_2=name0, beta_3=0, beta_4=0)
        print(f"pz_second(4,2,{name0},0,0,{name0},0,0): {fmt(v)}")

        if N >= 2:
            # A-B second derivative (off-diagonal A)
            v = partial_z_second(alpha_1=4, beta_1=2, basis_function=basis,
                                 alpha_2=name0, alpha_3=0, alpha_4=1,
                                 beta_2=name0, beta_3=0, beta_4=0)
            print(f"pz_second(4,2,{name0},0,1,{name0},0,0): {fmt(v)}")

        # Test calculate_C
        # C(u, u)
        v = calculate_C(alpha_1=1, alpha_2=name0, alpha_3=0, alpha_4=0,
                        beta_1=1, beta_2=name0, beta_3=0, beta_4=0,
                        basis_function=basis)
        print(f"C(1,{name0},0,0,1,{name0},0,0): {fmt(v)}")

        # C(B, B)
        v = calculate_C(alpha_1=2, alpha_2=name0, alpha_3=0, alpha_4=0,
                        beta_1=2, beta_2=name0, beta_3=0, beta_4=0,
                        basis_function=basis)
        print(f"C(2,{name0},0,0,2,{name0},0,0): {fmt(v)}")

        # C(R, B)
        v = calculate_C(alpha_1=3, alpha_2=name0, alpha_3=0, alpha_4=0,
                        beta_1=2, beta_2=name0, beta_3=0, beta_4=0,
                        basis_function=basis)
        print(f"C(3,{name0},0,0,2,{name0},0,0): {fmt(v)}")

        # C(A, B) diagonal A
        v = calculate_C(alpha_1=4, alpha_2=name0, alpha_3=0, alpha_4=0,
                        beta_1=2, beta_2=name0, beta_3=0, beta_4=0,
                        basis_function=basis)
        print(f"C(4,{name0},0,0,2,{name0},0,0): {fmt(v)}")


if __name__ == '__main__':
    main()
