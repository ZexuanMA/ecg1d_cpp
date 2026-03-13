#!/usr/bin/env python3
"""Validate Phase 3 (Hamiltonian gradients + TDVP) against Python reference."""
import sys, os, subprocess, tempfile
from pathlib import Path
import numpy as np

# Add sibling Python ECG directory to path (robust to custom HOME settings)
ecg_py_dir = Path(__file__).resolve().parent.parent / 'ecg1d'
sys.path.insert(0, str(ecg_py_dir))

from Parameters_z import BasisParams
from Calculate_Hamiltonian_Kinetic_partial import calculate_Hamiltonian_Kinetic_partial
from Calculate_Hamiltonian_Harmonic_partial import calculate_Hamiltonian_Harmonic_partial
from overlap import overlap
from Hamiltonian_Kinetic_energy_functional import kinetic_energy_functional
from Hamiltonian_Harmonic_functional import Harmonic_functional

dtype = np.complex128

# Same initial conditions as the C++ test (1-particle)
basis = []
basis.append(BasisParams.from_arrays(
    u=1.0,
    A=np.array([[0.1+0.1j]]).astype(dtype),
    B=np.array([[0.3+0.1j]]).astype(dtype),
    R=np.array([0.2+0.1j]).astype(dtype),
    dtype=dtype, name=0, enforce_shapes=True
))

N = 1
basis_n = len(basis)

print("=== Python Phase 3: Hamiltonian Gradients (1-particle) ===")

# Test gradient values
for a1 in [1, 2, 3]:
    for real_flag in [True, False]:
        gT = calculate_Hamiltonian_Kinetic_partial(a1, 0, 0, 0, real_flag, basis)
        gW = calculate_Hamiltonian_Harmonic_partial(a1, 0, 0, 0, real_flag, basis)
        g_total = gT + gW
        label = f"grad(a1={a1}, Real={real_flag})"
        print(f"  {label}: gT={gT:.15e}, gW={gW:.15e}, total={g_total:.15e}")

# Initial energy
S = overlap(basis)
T = kinetic_energy_functional(basis)
V = Harmonic_functional(basis)
E = (T + V) / S
print(f"\nInitial E = {E.real:.15e}")
print(f"Expected ground state = 0.5")

# Run TDVP for a few steps using same code as Verification_One_Particle_Harmonic_1.py
from Calculate_C import calculate_C
import copy

alpha_z_list = []
for i in range(basis_n):
    alpha_z_list.append((1, i, 0, 0))
for i in range(basis_n):
    for j in range(N):
        alpha_z_list.append((2, i, j, 0))
for i in range(basis_n):
    for j in range(N):
        alpha_z_list.append((3, i, j, 0))

def assemble_C_py(alpha_z_list, basis_function, *, reg=0):
    d = len(alpha_z_list)
    C_1 = np.zeros((d, d), dtype=np.complex128)
    for a, alpha in enumerate(alpha_z_list):
        a1,a2,a3,a4 = alpha
        for b, beta in enumerate(alpha_z_list):
            b1,b2,b3,b4 = beta
            C_1[a, b] = calculate_C(a1,a2,a3,a4, b1,b2,b3,b4, basis_function)
    return C_1

def grad_H_py(alpha, basis_function, Real):
    a1,a2,a3,a4 = alpha
    gT = calculate_Hamiltonian_Kinetic_partial(a1,a2,a3,a4, Real, basis_function)
    gW = calculate_Hamiltonian_Harmonic_partial(a1,a2,a3,a4, Real, basis_function)
    return gT + gW

def assemble_grad_py(alpha_z_list, basis_function, Real):
    g = np.zeros(len(alpha_z_list), dtype=np.complex128)
    for a, alpha in enumerate(alpha_z_list):
        g[a] = grad_H_py(alpha, basis_function, Real=Real)
    return g

# First TDVP step comparison
updata_constant = sum(1 for a in alpha_z_list if a[0] == 1)

C_1 = assemble_C_py(alpha_z_list, basis)
C_bar = C_1.conjugate().T

g_bar = assemble_grad_py(alpha_z_list, basis, Real=False)

C_bar_update = C_bar[updata_constant:, updata_constant:]
g_bar_update = g_bar[updata_constant:]

dz = np.linalg.lstsq(C_bar_update, -g_bar_update, rcond=1e-4)[0]

print(f"\nPython C_bar_update:")
print(C_bar_update)
print(f"\nPython g_bar_update: {g_bar_update}")
print(f"Python dz: {dz}")
print(f"|dz| = {np.linalg.norm(dz):.15e}")

# After 1 step
basis_trial = copy.deepcopy(basis)
dtao = 1e-3
for i in range(len(alpha_z_list)):
    a1,a2,a3,a4 = alpha_z_list[i]
    if a1 == 2:
        basis_trial[a2].B[a3][a3] += dz[i-updata_constant] * dtao
    if a1 == 3:
        basis_trial[a2].R[a3] += dz[i-updata_constant] * dtao

S_new = overlap(basis_trial)
T_new = kinetic_energy_functional(basis_trial)
V_new = Harmonic_functional(basis_trial)
E_new = (T_new + V_new) / S_new
print(f"\nAfter 1 step: E = {E_new.real:.15e}")
print(f"dE = {(E_new - E).real:.15e}")
