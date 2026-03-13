#!/usr/bin/env python3
"""Compare C++ outputs against Python reference values."""
import subprocess
import sys
import os
import re

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ecg1d'))

from Parameters_z import BasisParams
from Generate_basis_function import generate_basis_function
from overlap import overlap
from Hamiltonian_Kinetic_energy_functional import kinetic_energy_functional
from Hamiltonian_Harmonic_functional import Harmonic_functional
from Hamiltonian_Delta_contact_interation_functional import Delta_contact_interation_functional
from Hamiltonian_Gaussian_interaction_functional import Gaussian_interaction_functional
from Hamiltonian_Kicking_term_functional import kicking_term_functional


def parse_cpp_output(text):
    """Parse C++ output lines like 'Overlap:   1.23 + 4.56i'"""
    results = {}
    for line in text.strip().split('\n'):
        if ':' not in line or '===' in line or 'N=' in line:
            continue
        key, val = line.split(':', 1)
        key = key.strip()
        val = val.strip()
        # Parse "real + imagi"
        m = re.match(r'([-\d.e+]+)\s*\+\s*([-\d.e+]+)i', val)
        if m:
            results[key] = complex(float(m.group(1)), float(m.group(2)))
    return results


def run_test(N, basis_n, cpp_exe):
    csv_file = f"basis_N{N}_K{basis_n}.csv"
    if not os.path.exists(csv_file):
        print(f"  SKIP: {csv_file} not found")
        return 0, 0

    # Run C++
    result = subprocess.run([cpp_exe, '--csv', csv_file, str(N), str(basis_n)],
                           capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  FAIL: C++ crashed: {result.stderr}")
        return 0, 1

    cpp_vals = parse_cpp_output(result.stdout)

    # Compute Python reference
    basis = generate_basis_function(basis_n, N)
    py_vals = {}
    py_vals['Overlap'] = overlap(basis)
    py_vals['Kinetic'] = kinetic_energy_functional(basis)
    py_vals['Harmonic'] = Harmonic_functional(basis)
    if N >= 2:
        py_vals['Delta'] = Delta_contact_interation_functional(basis)
        py_vals['Gaussian'] = Gaussian_interaction_functional(basis)
    py_vals['Kicking'] = kicking_term_functional(basis)

    passed = 0
    failed = 0
    for key in py_vals:
        if key not in cpp_vals:
            print(f"  MISSING: {key}")
            failed += 1
            continue
        py_v = py_vals[key]
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

    for N in [1, 2, 3]:
        for basis_n in [1, 2, 3]:
            print(f"\n=== N={N}, K={basis_n} ===")
            p, f = run_test(N, basis_n, cpp_exe)
            total_pass += p
            total_fail += f

    print(f"\n{'='*40}")
    print(f"TOTAL: {total_pass} PASS, {total_fail} FAIL")


if __name__ == '__main__':
    main()
