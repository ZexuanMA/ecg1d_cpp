# Wavefunction Dynamics in `main_verification.cpp`

This note reviews the physics case implemented by
`main_verification.cpp` and derives the wavefunction dynamics used in the
Step-4 figures.

The central point is:

> Step 4 is a sudden-quench calculation.  The initial state is the ground
> state of a harmonic trap plus a cosine lattice.  At `t = 0` the cosine term
> is switched off, and the state evolves under the harmonic oscillator only.

So the plotted dynamics is not a repeated-kick calculation.  It is the free
harmonic-oscillator motion of a non-HO initial wavepacket.

---

## 1. What the code is solving

The constants are defined in `src/physical_constants.hpp`:

```cpp
hbar  = 1
mass  = 1
omega = 1
k_L   = 0.5
kappa = 1
```

Therefore

$$
2 k_L = 1,
$$

and all formulas below are in dimensionless oscillator units:

$$
\hbar = m = \omega = 1.
$$

### Step 1: prepare the initial state

`step1_imag_time_ecg.cpp` builds the ground state of

$$
\hat H_{\rm full}
= \hat T + \hat V_{\rm ho} + \hat V_{\rm cos}
= -\frac{1}{2}\frac{d^2}{dx^2}
  + \frac{1}{2}x^2
  + \cos x .
$$

The initial state is

$$
\psi_0(x)
= \psi(x,0)
\approx \text{ground state of } \hat H_{\rm full}.
$$

Because both \(x^2\) and \(\cos x\) are even functions, the exact ground state
is real and even:

$$
\psi_0(-x) = \psi_0(x).
$$

Small imaginary or odd components in the ECG output are numerical/variational
artifacts, not physical features of the exact Hamiltonian.

### Step 4: real-time dynamics after the quench

`step4_realtime_tdvp.cpp` uses

```cpp
HamiltonianTerms terms = HamiltonianTerms::kinetic_harmonic();
```

so the real-time Hamiltonian is

$$
\hat H_0
= -\frac{1}{2}\frac{d^2}{dx^2}
  + \frac{1}{2}x^2 .
$$

The exact Schrodinger equation after the quench is

$$
i\frac{\partial}{\partial t}\psi(x,t)
= \hat H_0 \psi(x,t),
\qquad
\psi(x,0) = \psi_0(x).
$$

The default total time in `main_verification.cpp` is

$$
T_{\rm total}=2\pi,
$$

which is one harmonic-oscillator period.

---

## 2. Solve the Schrodinger equation step by step

We solve

$$
i\partial_t\psi(x,t)
=
\left(
-\frac{1}{2}\partial_x^2 + \frac{1}{2}x^2
\right)\psi(x,t).
$$

### Step 2.1: solve the stationary harmonic oscillator

Look for separated solutions

$$
\psi_n(x,t) = \phi_n(x)e^{-iE_nt}.
$$

Substitute this into the Schrodinger equation:

$$
i\partial_t\left[\phi_n(x)e^{-iE_nt}\right]
=
\hat H_0\phi_n(x)e^{-iE_nt}.
$$

The left side is

$$
i(-iE_n)\phi_n(x)e^{-iE_nt}
=
E_n\phi_n(x)e^{-iE_nt}.
$$

Cancel the common time factor:

$$
\hat H_0\phi_n(x) = E_n\phi_n(x).
$$

This is the time-independent harmonic-oscillator eigenproblem:

$$
\left(
-\frac{1}{2}\frac{d^2}{dx^2}
+\frac{1}{2}x^2
\right)\phi_n(x)
= E_n\phi_n(x).
$$

Its normalized eigenfunctions are

$$
\phi_n(x)
=
\frac{1}{\pi^{1/4}\sqrt{2^n n!}}
H_n(x)e^{-x^2/2},
$$

where \(H_n(x)\) is the Hermite polynomial.  The eigenvalues are

$$
E_n = n + \frac{1}{2},
\qquad n=0,1,2,\ldots
$$

These \(\phi_n\) form a complete orthonormal basis.

### Step 2.2: expand the initial wavefunction

The initial state is not the harmonic-oscillator ground state.  It is the
ground state of

$$
\hat H_{\rm full} = \hat H_0 + \cos x.
$$

Therefore expand it in the harmonic-oscillator eigenbasis:

$$
\psi_0(x)
=
\sum_{n=0}^{\infty} c_n \phi_n(x),
$$

with coefficients

$$
c_n
=
\langle \phi_n|\psi_0\rangle
=
\int_{-\infty}^{\infty}\phi_n^*(x)\psi_0(x)\,dx.
$$

Since the exact \(\psi_0\) is even and \(\phi_n\) has parity

$$
\phi_n(-x)=(-1)^n\phi_n(x),
$$

all odd coefficients vanish:

$$
c_1=c_3=c_5=\cdots=0.
$$

Only even oscillator levels participate:

$$
\psi_0(x)
=
\sum_{m=0}^{\infty} c_{2m}\phi_{2m}(x).
$$

This fact is important for the figures: the wavepacket mainly shows a
breathing motion, not center-of-mass motion.

### Step 2.3: evolve each eigenmode

Because the Schrodinger equation is linear, every harmonic-oscillator
eigenmode evolves independently:

$$
\phi_n(x)
\longrightarrow
e^{-iE_nt}\phi_n(x).
$$

Therefore the exact solution is

$$
\boxed{
\psi(x,t)
=
\sum_{n=0}^{\infty}
c_n e^{-i(n+1/2)t}\phi_n(x)
}
$$

or, using even parity,

$$
\boxed{
\psi(x,t)
=
e^{-it/2}
\sum_{m=0}^{\infty}
c_{2m} e^{-i2mt}\phi_{2m}(x).
}
$$

The global phase \(e^{-it/2}\) does not affect density:

$$
n(x,t)=|\psi(x,t)|^2.
$$

The relative phases \(e^{-i2mt}\) are what change the shape of the packet.

### Step 2.4: periodicity

For a general harmonic-oscillator wavefunction, the full wavefunction returns
after \(2\pi\), up to a global phase:

$$
\psi(x,t+2\pi) = -\psi(x,t).
$$

The density is exactly periodic:

$$
n(x,t+2\pi) = n(x,t).
$$

For an exactly even initial state, only \(n=2m\) terms appear, so the relative
phase already returns after \(\pi\):

$$
e^{-i2m(t+\pi)} = e^{-i2mt}.
$$

Thus

$$
\psi(x,t+\pi) = e^{-i\pi/2}\psi(x,t),
$$

and therefore

$$
n(x,t+\pi)=n(x,t).
$$

So the exact density breathing period is \(\pi\), while the complete
oscillator period is \(2\pi\).  If the ECG initial state has a small odd or
complex component, the numerical trace can also show a tiny center-of-mass
oscillation with period \(2\pi\).

---

## 3. Coordinate-space propagator form

The same solution can be written directly in coordinate space:

$$
\psi(x,t)
=
\int_{-\infty}^{\infty}
K(x,x';t)\psi_0(x')\,dx',
$$

where the harmonic-oscillator propagator is

$$
K(x,x';t)
=
\frac{1}{\sqrt{2\pi i\sin t}}
\exp\left[
\frac{i}{2\sin t}
\left((x^2+x'^2)\cos t - 2xx'\right)
\right],
$$

for \(t\neq n\pi\).  At \(t=n\pi\), this kernel is understood by taking the
limit.

Special times:

- \(t=0\): \(\psi(x,0)=\psi_0(x)\).
- \(t=\pi/2\): the oscillator evolution is essentially a Fourier transform.
  Position and momentum distributions are exchanged.
- \(t=\pi\): \(\psi(x,\pi)=e^{-i\pi/2}\psi_0(-x)\).  For an even initial
  state this is the same shape.
- \(t=2\pi\): \(\psi(x,2\pi)=-\psi_0(x)\).  The density is exactly back to the
  initial density.

This explains why the position-density heatmap should breathe and then return,
rather than diffuse away.

---

## 4. Moment dynamics: the simplest way to read the figures

For the harmonic oscillator,

$$
\hat x(t)=\hat x(0)\cos t+\hat p(0)\sin t,
$$

$$
\hat p(t)=\hat p(0)\cos t-\hat x(0)\sin t.
$$

Therefore

$$
\langle x\rangle(t)
=
\langle x\rangle_0\cos t+\langle p\rangle_0\sin t,
$$

$$
\langle p\rangle(t)
=
\langle p\rangle_0\cos t-\langle x\rangle_0\sin t.
$$

For the exact even real initial ground state,

$$
\langle x\rangle_0 = 0,
\qquad
\langle p\rangle_0 = 0,
$$

so

$$
\langle x\rangle(t)=\langle p\rangle(t)=0.
$$

Any visible nonzero \(\langle x\rangle\) or \(\langle p\rangle\) in the ECG
trace is therefore a numerical/variational asymmetry, not a property of the
exact quench.

For the width,

$$
\langle x^2\rangle(t)
=
\langle x^2\rangle_0\cos^2 t
+\langle p^2\rangle_0\sin^2 t
+\frac{1}{2}\langle xp+px\rangle_0\sin 2t.
$$

Similarly,

$$
\langle p^2\rangle(t)
=
\langle p^2\rangle_0\cos^2 t
+\langle x^2\rangle_0\sin^2 t
-\frac{1}{2}\langle xp+px\rangle_0\sin 2t.
$$

For an even real initial state,

$$
\langle xp+px\rangle_0=0,
$$

so

$$
\boxed{
\langle x^2\rangle(t)
=
\langle x^2\rangle_0\cos^2 t
+\langle p^2\rangle_0\sin^2 t
}
$$

and

$$
\boxed{
\langle p^2\rangle(t)
=
\langle p^2\rangle_0\cos^2 t
+\langle x^2\rangle_0\sin^2 t.
}
$$

Thus the position width and momentum width exchange every quarter period.

From the current `out/verify/step1_ecg_energy_N1_K5.csv`:

```text
E_kinetic  = 0.145331047104415
E_harmonic = 0.437862782808295
```

Since

$$
E_{\rm kin}=\frac{1}{2}\langle p^2\rangle,
\qquad
E_{\rm ho}=\frac{1}{2}\langle x^2\rangle,
$$

the initial moments are approximately

$$
\langle p^2\rangle_0 \approx 0.290662,
\qquad
\langle x^2\rangle_0 \approx 0.875726.
$$

The initial packet is much wider in position than the harmonic-oscillator
ground state, for which

$$
\langle x^2\rangle = \langle p^2\rangle = \frac{1}{2}.
$$

Therefore the expected motion is:

- at \(t=0\): broad position density, narrow momentum density;
- at \(t=\pi/2\): narrow position density, broad momentum density;
- at \(t=\pi\): broad position density again;
- at \(t=3\pi/2\): narrow position density again;
- at \(t=2\pi\): return to the initial density.

This is the main physical interpretation of the Step-4 density and momentum
figures.

The conserved post-quench energy is

$$
E_0
=
\frac{1}{2}\left(\langle p^2\rangle+\langle x^2\rangle\right)
\approx
\frac{1}{2}(0.290662+0.875726)
\approx 0.583194.
$$

That matches the Step-4 ECG snapshot energy:

```text
out/verify/step4_ecg_snap_N1_K5.csv:
E(t) = 0.5831938299...
```

This energy is lower than the Step-1 full-Hamiltonian energy

$$
E_{\rm full}
\approx 1.22264,
$$

because the cosine contribution is removed at the quench:

$$
E_{\rm full}
= E_{\rm kin}+E_{\rm ho}+E_{\rm cos},
\qquad
E_0
= E_{\rm kin}+E_{\rm ho}.
$$

---

## 5. Momentum-space wavefunction

The code uses the Fourier convention

$$
\tilde\psi(k,t)
=
\frac{1}{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\psi(x,t)e^{-ikx}\,dx.
$$

The momentum distribution plotted as `n_k` is

$$
n(k,t)=|\tilde\psi(k,t)|^2.
$$

In the harmonic oscillator, the Fourier transform of an eigenfunction is the
same eigenfunction up to a phase:

$$
\mathcal{F}[\phi_n](k)=(-i)^n\phi_n(k).
$$

Therefore

$$
\tilde\psi(k,t)
=
\sum_{n=0}^{\infty}
c_n e^{-i(n+1/2)t}(-i)^n\phi_n(k).
$$

This is the momentum-space version of the same mode beating.  At a quarter
period, the position and momentum distributions are swapped.  That is why the
position-density heatmap and the momentum-density heatmap should look like
phase-space rotations of each other.

---

## 6. How the reference solver implements the exact solution

`step4_reference_dynamics.cpp` implements the finite-dimensional version of
the spectral derivation above.

### Step 6.1: build the post-quench Hamiltonian

The code builds

$$
H_{\rm evolve}
=
T_{\rm grid/DVR} + \operatorname{diag}\left(\frac{1}{2}x_j^2\right).
$$

The cosine term is not included.

### Step 6.2: diagonalize once

The matrix eigenproblem is

$$
H_{\rm evolve} U = U\Lambda,
$$

where

$$
\Lambda = \operatorname{diag}(\lambda_0,\lambda_1,\ldots).
$$

In code:

```cpp
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
lambda = es.eigenvalues();
U      = es.eigenvectors();
```

### Step 6.3: project the ECG initial wavefunction onto the reference basis

The ECG initial state is evaluated on the grid and normalized:

```cpp
psi0 = ecg_wavefunction_1p(ecg_psi0_basis, x);
normalize_grid_state(psi0, dx);
```

Then

$$
c = U^T\psi_0.
$$

### Step 6.4: propagate exactly in that basis

For every requested time,

$$
\psi(t)
=
U \operatorname{diag}\left(e^{-i\lambda_n t}\right)c.
$$

This is exactly the same formula as

$$
\psi(x,t)=\sum_n c_n e^{-iE_nt}\phi_n(x),
$$

but using the numerical grid/DVR eigenfunctions instead of analytic Hermite
functions.

The reference solver then writes:

- `step4_grid_density_N1.csv`: \(n(x,t)=|\psi(x,t)|^2\);
- `step4_grid_nk_N1.csv`: \(n(k,t)=|\tilde\psi(k,t)|^2\);
- `step4_grid_trace_N1.csv`: \(E(t)\), norm, \(\langle x\rangle\),
  \(\langle p\rangle\), fidelity.

The DVR path is the same, with the sinc-DVR kinetic matrix.

---

## 7. How the ECG TDVP solver approximates the same equation

For \(N=1\), the ECG wavefunction used by
`src/verify/ecg_wavefunction.cpp` is

$$
\psi_{\rm ECG}(x,t)
=
\sum_{i=1}^{K}
u_i(t)
\exp\left[
-(A_i(t)+B_i(t))x^2
+2R_i(t)B_i(t)x
-R_i(t)B_i(t)R_i(t)
\right].
$$

Collect all variational parameters into

$$
z=(u,A,B,R).
$$

The exact Schrodinger velocity is

$$
|\dot\psi\rangle_{\rm exact}
=
-i\hat H_0|\psi\rangle.
$$

But the ECG ansatz can only move inside its variational manifold.  TDVP chooses
the closest velocity available in the tangent space:

$$
|\dot\psi\rangle_{\rm TDVP}
=
\sum_\beta \dot z_\beta |\partial_\beta\psi\rangle.
$$

Projecting the Schrodinger equation onto the tangent vectors gives

$$
i\sum_\beta
C_{\alpha\beta}\dot z_\beta
=
g_\alpha,
$$

where, schematically,

$$
C_{\alpha\beta}
\sim
\langle \partial_\alpha\psi|
\partial_\beta\psi\rangle,
\qquad
g_\alpha
\sim
\langle \partial_\alpha\psi|\hat H_0|\psi\rangle.
$$

So

$$
\boxed{
\dot z = -i C^{-1}g.
}
$$

That is exactly the sign structure documented in
`src/realtime_tdvp.hpp`:

```cpp
// real time: i C z_dot = g -> dz = -i C^{-1} g
```

The implementation uses RK4 and an SVD pseudoinverse/regularization of \(C\).
This is a variational approximation to the exact spectral solution, not a new
physical Hamiltonian.

### Important diagnostic

The exact harmonic oscillator conserves norm and energy:

$$
\frac{d}{dt}\langle\psi|\psi\rangle=0,
\qquad
\frac{d}{dt}\langle H_0\rangle=0.
$$

The reference grid/DVR traces should therefore have norm \(1\) and constant
energy up to roundoff.

If the raw ECG TDVP trace shows norm drift, that is not physical.  It comes
from the projected/regularized TDVP dynamics, for example SVD truncation or
Tikhonov regularization in the tangent-space solve.  The current
`step4_ecg_trace_N1_K5.csv` shows this kind of norm drift, while the normalized
energy remains nearly constant.  This is a TDVP numerical issue, not part of
the exact Schrodinger solution.

---

## 8. Interpretation of the physics figures

### Density heatmap \(n(x,t)\)

This should show a breathing packet:

1. The initial cosine-on ground state is broad in position.
2. After the cosine is removed, the harmonic trap rotates the wavepacket in
   phase space.
3. At \(t=\pi/2\), the position distribution becomes narrow.
4. At \(t=\pi\), the position distribution returns to its initial broad shape.
5. The same pattern repeats until \(2\pi\).

There should be no physical drift of the packet center for the exact even
initial state.

### Momentum heatmap \(n(k,t)\)

The momentum distribution breathes out of phase with the position density:

1. At \(t=0\), the packet is broad in position and narrow in momentum.
2. At \(t=\pi/2\), the packet is narrow in position and broad in momentum.
3. At \(t=\pi\), the initial momentum distribution returns.

This follows directly from the oscillator phase-space rotation.

### Energy plot

After the quench, the conserved energy is

$$
E_0
=
\langle\psi_0|\hat H_0|\psi_0\rangle,
$$

not

$$
\langle\psi_0|\hat H_{\rm full}|\psi_0\rangle.
$$

So the Step-4 energy should be around

$$
E_0 \approx 0.583,
$$

while the Step-1 full energy is around

$$
E_{\rm full}\approx 1.223.
$$

The difference is the removed cosine potential contribution.

### Fidelity plot

The fidelity with the initial state is

$$
F(t)
=
\left|
\langle\psi_0|\psi(t)\rangle
\right|^2.
$$

Using the HO expansion,

$$
F(t)
=
\left|
\sum_n |c_n|^2 e^{-i(n+1/2)t}
\right|^2.
$$

For an exactly even state,

$$
F(t)
=
\left|
\sum_m |c_{2m}|^2 e^{-i(2m+1/2)t}
\right|^2.
$$

The global phase cancels in the absolute square.  Therefore \(F(t)\) is also
periodic with period \(\pi\) for a purely even initial state.  It is smallest
when the relative phases of the occupied even HO modes are most dephased.

### Wavefunction snapshots

The complex wavefunction evolves by phase interference between HO eigenmodes.
The density can return while the wavefunction has acquired a global phase.
For example:

$$
\psi(x,\pi)=e^{-i\pi/2}\psi_0(x)
$$

for an exact even initial state, and

$$
\psi(x,2\pi)=-\psi_0(x).
$$

The sign or global phase does not change density.

---

## 9. N-particle extension for the non-interacting case

The Step-4 reference solver in `main_verification.cpp` is currently written
only for `N = 1`.  However, the Hamiltonian used here has no inter-particle
interaction:

$$
\hat H_N
=
\sum_{a=1}^{N}
\left(
-\frac{1}{2}\frac{\partial^2}{\partial x_a^2}
+\frac{1}{2}x_a^2
\right).
$$

If the initial many-body state is a bosonic product state,

$$
\Psi(x_1,\ldots,x_N,0)
=
\prod_{a=1}^{N}\psi_0(x_a),
$$

then the time-dependent state stays factorized:

$$
\boxed{
\Psi(x_1,\ldots,x_N,t)
=
\prod_{a=1}^{N}\psi(x_a,t).
}
$$

The one-body density is still

$$
n(x,t)=|\psi(x,t)|^2,
$$

and the total energy is

$$
E_N = N E_1.
$$

This is why the `N=2` non-interacting verification can be understood by first
solving the single-particle problem.

---

## 10. Practical checklist when reading the current outputs

1. If the figure is from `step4_grid_*` or `step4_dvr_*`, read it as the exact
   post-quench harmonic-oscillator dynamics.
2. If the figure is from `step4_ecg_*`, remember it is the TDVP approximation
   to the same dynamics.
3. The physical energy after the quench is \(E_{\rm kin}+E_{\rm ho}\), not the
   Step-1 full energy.
4. The density should breathe with period \(\pi\) for an exact even initial
   state.
5. The full wavefunction has oscillator phase factors, so a returned density
   does not necessarily mean the complex wavefunction is identical without a
   global phase.
6. Norm drift in ECG TDVP is numerical.  It is not predicted by the exact
   Schrodinger equation.

### Code-review caveat: wavefunction fidelity from density is not exact

`step4_cross_check.cpp` compares densities correctly, but its approximate
N=1 reference wavefunction fidelity reconstructs the reference wavefunction as

$$
\psi_{\rm ref}(x,t) \approx \sqrt{n_{\rm ref}(x,t)}
$$

with zero phase.  That is not the exact real-time wavefunction in general.
Even if the Hamiltonian matrix and initial state are real, the unitary factor
\(e^{-iHt}\) creates a complex phase pattern unless the initial state is a
single eigenstate or the time is a special recurrence time.

Therefore:

- density metrics such as \(L^2[n_x]\), \(L^2[n_k]\), and Bhattacharyya
  overlap are the reliable Step-4 cross-check diagnostics;
- the complex wavefunction fidelity is reliable only if the reference complex
  \(\psi(x,t)\) is saved from the spectral propagation, or if the time is a
  special point where the wavefunction is known to be real up to a global
  phase;
- this caveat does not affect the exact derivation above, only the way one
  diagnostic is reconstructed from saved output.
