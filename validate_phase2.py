#!/usr/bin/env python3
"""Compare C++ Phase 2 (derivatives) outputs against Python reference values."""
import subprocess
import sys
import os
import re

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'ecg1d'))

from Parameters_z import BasisParams
from Generate_basis_function import generate_basis_function
from Partial_z_first import partial_z_first
from Partial_z_second import partial_z_second
from Calculate_C import calculate_C
from overlap import overlap
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
    name0 = basis[0].name
    name1 = basis[1].name

    py_vals = {}

    # Phase 1 (quick sanity check)
    py_vals['Overlap'] = overlap(basis)

    # partial_z_first
    py_vals[f'pz_first(1,T,{name0},0,0)'] = partial_z_first(1, True, basis, alpha_2=name0, alpha_3=0, alpha_4=0)
    py_vals[f'pz_first(1,F,{name0},0,0)'] = partial_z_first(1, False, basis, alpha_2=name0, alpha_3=0, alpha_4=0)
    py_vals[f'pz_first(2,T,{name0},0,0)'] = partial_z_first(2, True, basis, alpha_2=name0, alpha_3=0, alpha_4=0)
    py_vals[f'pz_first(2,F,{name0},0,0)'] = partial_z_first(2, False, basis, alpha_2=name0, alpha_3=0, alpha_4=0)
    py_vals[f'pz_first(3,T,{name0},0,0)'] = partial_z_first(3, True, basis, alpha_2=name0, alpha_3=0, alpha_4=0)
    py_vals[f'pz_first(3,F,{name1},0,0)'] = partial_z_first(3, False, basis, alpha_2=name1, alpha_3=0, alpha_4=0)
    py_vals[f'pz_first(4,T,{name0},0,0)'] = partial_z_first(4, True, basis, alpha_2=name0, alpha_3=0, alpha_4=0)
    if N >= 2:
        py_vals[f'pz_first(4,T,{name0},0,1)'] = partial_z_first(4, True, basis, alpha_2=name0, alpha_3=0, alpha_4=1)

    # partial_z_second
    py_vals[f'pz_second(2,2,{name0},0,0,{name0},0,0)'] = partial_z_second(2, 2, basis, alpha_2=name0, alpha_3=0, alpha_4=0, beta_2=name0, beta_3=0, beta_4=0)
    py_vals[f'pz_second(1,1,{name0},0,0,{name0},0,0)'] = partial_z_second(1, 1, basis, alpha_2=name0, alpha_3=0, alpha_4=0, beta_2=name0, beta_3=0, beta_4=0)
    py_vals[f'pz_second(1,2,{name0},0,0,{name0},0,0)'] = partial_z_second(1, 2, basis, alpha_2=name0, alpha_3=0, alpha_4=0, beta_2=name0, beta_3=0, beta_4=0)
    py_vals[f'pz_second(3,2,{name0},0,0,{name0},0,0)'] = partial_z_second(3, 2, basis, alpha_2=name0, alpha_3=0, alpha_4=0, beta_2=name0, beta_3=0, beta_4=0)
    py_vals[f'pz_second(4,2,{name0},0,0,{name0},0,0)'] = partial_z_second(4, 2, basis, alpha_2=name0, alpha_3=0, alpha_4=0, beta_2=name0, beta_3=0, beta_4=0)
    if N >= 2:
        py_vals[f'pz_second(4,2,{name0},0,1,{name0},0,0)'] = partial_z_second(4, 2, basis, alpha_2=name0, alpha_3=0, alpha_4=1, beta_2=name0, beta_3=0, beta_4=0)

    # calculate_C
    py_vals[f'C(1,{name0},0,0,1,{name0},0,0)'] = calculate_C(1, name0, 0, 0, 1, name0, 0, 0, basis)
    py_vals[f'C(2,{name0},0,0,2,{name0},0,0)'] = calculate_C(2, name0, 0, 0, 2, name0, 0, 0, basis)
    py_vals[f'C(3,{name0},0,0,2,{name0},0,0)'] = calculate_C(3, name0, 0, 0, 2, name0, 0, 0, basis)
    py_vals[f'C(4,{name0},0,0,2,{name0},0,0)'] = calculate_C(4, name0, 0, 0, 2, name0, 0, 0, basis)

    passed = 0
    failed = 0
    for key in py_vals:
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
        basis_n = 2
        print(f"\n=== N={N}, K={basis_n} ===")
        p, f = run_test(N, basis_n, cpp_exe)
        total_pass += p
        total_fail += f

    print(f"\n{'='*40}")
    print(f"TOTAL: {total_pass} PASS, {total_fail} FAIL")


if __name__ == '__main__':
    main()
