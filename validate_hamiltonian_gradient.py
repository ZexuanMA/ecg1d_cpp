#!/usr/bin/env python3
"""Compare C++ Hamiltonian gradient (first derivative) outputs against Python reference values."""
import subprocess
import sys
import os
import re

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ecg1d'))

from Parameters_z import BasisParams
from Generate_basis_function import generate_basis_function
from Calculate_Hamiltonian_Kinetic_partial import calculate_Hamiltonian_Kinetic_partial
from Calculate_Hamiltonian_Harmonic_partial import calculate_Hamiltonian_Harmonic_partial
from Calculate_Hamiltonian_Delta_partial import calculate_Hamiltonian_Delta_partial
from Calculate_Hamiltonian_Gaussian_partial import calculate_Hamiltonian_Gaussian_partial
from Calculate_Hamiltonian_Kicking_partial import calculate_Hamiltonian_Kicking_partial
import numpy as np


def parse_cpp_output(text):
    """Parse C++ output lines like 'label: 1.23 + 4.56i'"""
    results = {}
    for line in text.strip().split('\n'):
        if ':' not in line or '===' in line or 'N=' in line:
            continue
        key, val = line.split(':', 1)
        key = key.strip()
        val = val.strip()
        m = re.match(r'([-\d.e+]+)\s*\+\s*([-\d.e+]+)i', val)
        if m:
            results[key] = complex(float(m.group(1)), float(m.group(2)))
    return results


# Map from label prefix to Python function
GRAD_FUNCS = {
    'gKinetic':  calculate_Hamiltonian_Kinetic_partial,
    'gHarmonic': calculate_Hamiltonian_Harmonic_partial,
    'gDelta':    calculate_Hamiltonian_Delta_partial,
    'gGaussian': calculate_Hamiltonian_Gaussian_partial,
    'gKicking':  calculate_Hamiltonian_Kicking_partial,
}


def build_py_vals(basis, N):
    """Compute all Hamiltonian gradient reference values matching C++ output labels."""
    py_vals = {}

    name0 = basis[0].name
    name1 = basis[1].name if len(basis) > 1 else name0

    def rf_str(r):
        return 'T' if r else 'F'

    # Determine which terms to test
    terms = ['gKinetic', 'gHarmonic', 'gKicking']
    if N >= 2:
        terms += ['gDelta', 'gGaussian']

    # a1={1,2,3}, Real={T,F}, alpha_2=name0
    for a1 in [1, 2, 3]:
        for real_flag in [True, False]:
            tag = f"({a1},{rf_str(real_flag)},{name0},0,0)"
            for term in terms:
                fn = GRAD_FUNCS[term]
                py_vals[f"{term}{tag}"] = fn(a1, name0, 0, 0, real_flag, basis)

    # a1={2,3} with alpha_2=name1
    if len(basis) > 1:
        for a1 in [2, 3]:
            for real_flag in [True, False]:
                tag = f"({a1},{rf_str(real_flag)},{name1},0,0)"
                for term in terms:
                    fn = GRAD_FUNCS[term]
                    py_vals[f"{term}{tag}"] = fn(a1, name1, 0, 0, real_flag, basis)

    # a1=4 (A matrix derivative) if N >= 2
    if N >= 2:
        for real_flag in [True, False]:
            for a3, a4 in [(0, 0), (0, 1)]:
                tag = f"(4,{rf_str(real_flag)},{name0},{a3},{a4})"
                for term in terms:
                    fn = GRAD_FUNCS[term]
                    py_vals[f"{term}{tag}"] = fn(4, name0, a3, a4, real_flag, basis)

    return py_vals


def run_test(N, basis_n, cpp_exe):
    csv_file = f"basis_N{N}_K{basis_n}.csv"
    if not os.path.exists(csv_file):
        print(f"  SKIP: {csv_file} not found")
        return 0, 0

    result = subprocess.run([cpp_exe, '--csv', csv_file, str(N), str(basis_n)],
                           capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  FAIL: C++ crashed: {result.stderr}")
        return 0, 1

    cpp_vals = parse_cpp_output(result.stdout)
    basis = generate_basis_function(basis_n, N)
    py_vals = build_py_vals(basis, N)

    passed = 0
    failed = 0
    for key in sorted(py_vals):
        if key not in cpp_vals:
            print(f"  MISSING: {key}")
            failed += 1
            continue
        py_v = complex(py_vals[key])
        cpp_v = cpp_vals[key]
        rel_err = abs(cpp_v - py_v) / max(abs(py_v), 1e-15)
        status = "PASS" if rel_err < 1e-10 else "FAIL"
        if status == "FAIL":
            print(f"  {status}: {key}: py={py_v:.15e}, cpp={cpp_v:.15e}, rel_err={rel_err:.2e}")
            failed += 1
        else:
            print(f"  {status}: {key} (rel_err={rel_err:.2e})")
            passed += 1

    return passed, failed


def main():
    cpp_exe = './build/ecg1d'
    if not os.path.exists(cpp_exe):
        print("ERROR: Build C++ first: cd build && cmake .. && make")
        return

    total_pass = 0
    total_fail = 0

    for N in [1, 2]:
        for basis_n in [2, 3]:
            print(f"\n=== N={N}, K={basis_n} ===")
            p, f = run_test(N, basis_n, cpp_exe)
            total_pass += p
            total_fail += f

    print(f"\n{'='*40}")
    print(f"TOTAL: {total_pass} PASS, {total_fail} FAIL")


if __name__ == '__main__':
    main()
