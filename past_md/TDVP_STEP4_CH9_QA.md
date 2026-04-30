# Q&A: Understanding Chapter 9 of the Step-4 TDVP Report

This note follows the mathematical diagnosis in:

`out/verify/N=1_step4_verification.html`, Section 9.

The goal is to explain the formulas slowly, question by question.

---

## Question 1

Original question:

> How should I understand the tangent space in line 557?

Refined question:

> In Section 9.1, what does the tangent space
> 
> $$
> T_\psi\mathcal{M}
> =
> \mathrm{span}\{|\partial_\alpha\psi\rangle\}_{\alpha=1}^{d}
> $$
> 
> mean, and why is it important for TDVP?

## Answer

The ECG ansatz does not describe every possible wavefunction. It describes only
wavefunctions that can be written using the chosen ECG parameters.

The report writes this set as

$$
\mathcal{M}
=
\{\,|\psi(z)\rangle : z\in\mathbb{C}^{d}\,\}.
$$

Here

$$
z=(u,A,B,R)
$$

means all variational parameters of the ECG wavefunction. So
\(\mathcal{M}\) is the variational manifold: the set of all wavefunctions that
the ECG ansatz can represent.

Now suppose the current wavefunction is

$$
|\psi\rangle
=
|\psi(z)\rangle.
$$

If we make a very small change to the parameters,

$$
z
\longrightarrow
z+\delta z,
$$

then the wavefunction changes as

$$
|\psi(z+\delta z)\rangle
=
|\psi(z)\rangle
+
\sum_{\alpha=1}^{d}
\delta z_\alpha
|\partial_\alpha\psi\rangle
+
O(\|\delta z\|^2).
$$

The first-order change is therefore

$$
|\delta\psi\rangle
=
\sum_{\alpha=1}^{d}
\delta z_\alpha
|\partial_\alpha\psi\rangle.
$$

Each vector

$$
|\partial_\alpha\psi\rangle
=
\frac{\partial|\psi(z)\rangle}{\partial z_\alpha}
$$

is one possible infinitesimal direction in which the wavefunction can move when
we change one ECG parameter.

Therefore the tangent space is

$$
T_\psi\mathcal{M}
=
\left\{
\sum_{\alpha=1}^{d}
c_\alpha
|\partial_\alpha\psi\rangle
:
c_\alpha\in\mathbb{C}
\right\}.
$$

In words:

> The tangent space is the set of all infinitesimal wavefunction changes that
> can be produced by infinitesimally changing the ECG parameters at the current
> state.

For the ECG expansion

$$
|\psi\rangle
=
\sum_{k=1}^{K}
u_k|\phi_k(A,B,R)\rangle,
$$

one simple tangent direction comes from changing a linear coefficient \(u_k\):

$$
\frac{\partial|\psi\rangle}{\partial u_k}
=
|\phi_k\rangle.
$$

Other tangent directions come from changing the nonlinear Gaussian parameters:

$$
\frac{\partial|\psi\rangle}{\partial A},
\qquad
\frac{\partial|\psi\rangle}{\partial B},
\qquad
\frac{\partial|\psi\rangle}{\partial R}.
$$

These directions are more complicated, but the idea is the same: each one says
how the wavefunction changes when one parameter is nudged.

## Geometric Picture

Think of the full Hilbert space as a very large space. The ECG ansatz is only a
curved surface inside that space:

$$
\mathcal{M}
\subset
\mathcal{H}.
$$

At the current point \(|\psi\rangle\), the tangent space is like the tangent
plane touching that surface:

$$
T_\psi\mathcal{M}
\quad
\text{is the local linear approximation to}
\quad
\mathcal{M}.
$$

This tangent space is not the whole manifold. It is only the local set of
allowed velocities at the current wavefunction.

## Why TDVP Needs This

The exact Schrodinger equation wants the wavefunction to move with velocity

$$
|\dot\psi\rangle_{\mathrm{exact}}
=
-iH|\psi\rangle.
$$

But this exact velocity may point outside the ECG manifold:

$$
-iH|\psi\rangle
\notin
T_\psi\mathcal{M}.
$$

TDVP says: choose the closest velocity that the ECG ansatz can actually
produce. That means projecting the exact Schrodinger velocity onto the tangent
space:

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
-iP_T H|\psi\rangle.
$$

Here \(P_T\) is the orthogonal projector onto \(T_\psi\mathcal{M}\).

So the tangent space is the key object because TDVP does not evolve in the
whole Hilbert space. It evolves using only the directions available from the
current ECG parameters.

## Common Misunderstanding

The tangent space is not the set of all ECG wavefunctions. That set is
\(\mathcal{M}\).

The tangent space is also not fixed forever. It depends on the current
wavefunction:

$$
T_\psi\mathcal{M}
\quad
\text{changes when}
\quad
|\psi\rangle
\text{ changes}.
$$

The tangent space is only the first-order, local set of possible movements
available to the ansatz at one instant.

---

## Question 2

Original question:

> How does the Dirac-Frenkel principle minimise the residual? I need more
> mathematical details.

Refined question:

> In Section 9.1, why does minimising
> 
> $$
> \bigl\|(i\partial_t-H)|\psi\rangle\bigr\|^2
> $$
> 
> over tangent-space velocities lead to
> 
> $$
> i\,C_{\alpha\beta}\,\dot z_\beta
> =
> g_\alpha,
> \qquad
> C_{\alpha\beta}
> =
> \langle\partial_\alpha\psi|\partial_\beta\psi\rangle,
> \qquad
> g_\alpha
> =
> \langle\partial_\alpha\psi|H|\psi\rangle?
> $$

## Answer

Start from the exact Schrodinger equation:

$$
i\partial_t|\psi\rangle
=
H|\psi\rangle.
$$

Equivalently,

$$
(i\partial_t-H)|\psi\rangle
=
0.
$$

If the variational ansatz were exact, we could make this residual zero. But
inside a finite ECG manifold, the wavefunction can only move along the tangent
space:

$$
|\dot\psi\rangle
\in
T_\psi\mathcal{M}.
$$

Using the tangent basis,

$$
|\dot\psi\rangle
=
\sum_{\beta=1}^{d}
\dot z_\beta
|\partial_\beta\psi\rangle.
$$

The residual is therefore

$$
|R\rangle
=
i|\dot\psi\rangle
-
H|\psi\rangle
$$

or

$$
|R\rangle
=
i
\sum_{\beta=1}^{d}
\dot z_\beta
|\partial_\beta\psi\rangle
-
H|\psi\rangle.
$$

The Dirac-Frenkel principle chooses the velocity coefficients
\(\dot z_\beta\) that make this residual as small as possible:

$$
\min_{\dot z_1,\dots,\dot z_d}
\|R\|^2.
$$

This is a Hilbert-space least-squares problem.

## Matrix Form

Introduce the linear map \(T\) whose columns are the tangent vectors:

$$
T
=
\begin{bmatrix}
|\partial_1\psi\rangle
&
|\partial_2\psi\rangle
&
\cdots
&
|\partial_d\psi\rangle
\end{bmatrix}.
$$

Then

$$
|\dot\psi\rangle
=
T\dot z,
$$

and the residual becomes

$$
|R\rangle
=
iT\dot z
-
H|\psi\rangle.
$$

Define

$$
|b\rangle
=
H|\psi\rangle.
$$

Then the minimisation problem is

$$
\min_{\dot z}
\|iT\dot z-|b\rangle\|^2.
$$

This has exactly the same structure as an ordinary least-squares problem:

$$
\min_x
\|Ax-b\|^2.
$$

Here the role of \(A\) is played by \(iT\), and the role of \(x\) is played by
\(\dot z\).

## Expanding the Norm

The squared residual is

$$
\|R\|^2
=
\langle R|R\rangle.
$$

Substitute

$$
|R\rangle
=
iT\dot z
-
|b\rangle.
$$

Then

$$
\|R\|^2
=
\langle iT\dot z-b|iT\dot z-b\rangle.
$$

Using

$$
C
=
T^\dagger T,
\qquad
g
=
T^\dagger |b\rangle
=
T^\dagger H|\psi\rangle,
$$

we get

$$
\|R\|^2
=
\dot z^\dagger C\dot z
+
i\dot z^\dagger g
-
ig^\dagger\dot z
+
\langle b|b\rangle.
$$

The last term \(\langle b|b\rangle\) does not depend on \(\dot z\). To minimise
the residual, differentiate with respect to the complex conjugate variables
\(\dot z^\ast\):

$$
\frac{\partial}{\partial \dot z^\ast}
\|R\|^2
=
C\dot z
+
ig.
$$

Setting this derivative to zero gives

$$
C\dot z
=
-ig.
$$

Equivalently,

$$
iC\dot z
=
g.
$$

In index notation, this is exactly

$$
i\,C_{\alpha\beta}\,\dot z_\beta
=
g_\alpha.
$$

## Orthogonality Interpretation

There is an equivalent and more geometric way to say the same thing.

At the minimum, the residual must be orthogonal to every tangent vector:

$$
\langle\partial_\alpha\psi|R\rangle
=
0
\qquad
\text{for every } \alpha.
$$

Substitute the residual:

$$
\left\langle
\partial_\alpha\psi
\middle|
i\sum_{\beta}
\dot z_\beta
|\partial_\beta\psi\rangle
-
H|\psi\rangle
\right\rangle
=
0.
$$

Therefore

$$
i
\sum_{\beta}
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle
\dot z_\beta
-
\langle\partial_\alpha\psi|H|\psi\rangle
=
0.
$$

Using the definitions

$$
C_{\alpha\beta}
=
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle,
\qquad
g_\alpha
=
\langle\partial_\alpha\psi|H|\psi\rangle,
$$

we again obtain

$$
iC_{\alpha\beta}\dot z_\beta
=
g_\alpha.
$$

So the phrase "Dirac-Frenkel minimises the residual" means:

> Among all velocities that the ECG parameters can produce, choose the one for
> which the Schrodinger-equation error is orthogonal to the tangent space.

## Connection to Projection

If \(C\) is invertible, then

$$
\dot z
=
-iC^{-1}g.
$$

Since

$$
g
=
T^\dagger H|\psi\rangle,
$$

the wavefunction velocity is

$$
|\dot\psi\rangle
=
T\dot z
=
-iT C^{-1}T^\dagger H|\psi\rangle.
$$

Define

$$
P_T
=
TC^{-1}T^\dagger.
$$

Then

$$
|\dot\psi\rangle
=
-iP_T H|\psi\rangle.
$$

This is why TDVP is described as projected Schrodinger evolution. It takes the
exact velocity

$$
-iH|\psi\rangle
$$

and projects it onto the tangent space:

$$
-iH|\psi\rangle
\longrightarrow
-iP_T H|\psi\rangle.
$$

The part that cannot be represented by the current ECG tangent space is the
remaining residual:

$$
|R\rangle
=
P_T H|\psi\rangle
-
H|\psi\rangle
=
-(1-P_T)H|\psi\rangle.
$$

This residual is not zero in general, but it is orthogonal to the tangent space:

$$
T^\dagger |R\rangle
=
0.
$$

---

## Question 3

Original question:

> How do we make sure that the first-order change stays in the ECG space, or
> the variational manifold?

Refined question:

> When we write
> 
> $$
> |\psi(z+\delta z)\rangle
> =
> |\psi(z)\rangle
> +
> \sum_\alpha
> \delta z_\alpha
> |\partial_\alpha\psi\rangle
> +
> O(\|\delta z\|^2),
> $$
> 
> does the first-order change
> 
> $$
> \sum_\alpha
> \delta z_\alpha
> |\partial_\alpha\psi\rangle
> $$
> 
> itself stay inside the ECG variational manifold?

## Answer

The short answer is:

> The exact changed state \( |\psi(z+\delta z)\rangle \) stays inside the ECG
> manifold, but the first-order change by itself is a tangent vector, not
> usually an ECG wavefunction.

This distinction is very important.

The ECG manifold is

$$
\mathcal{M}
=
\{\,|\psi(z)\rangle : z\in\mathbb{C}^d\,\}.
$$

If we choose a new parameter vector

$$
z+\delta z,
$$

then by definition

$$
|\psi(z+\delta z)\rangle
\in
\mathcal{M}.
$$

So the full updated ansatz state is still an ECG state.

But the Taylor expansion says

$$
|\psi(z+\delta z)\rangle
=
|\psi(z)\rangle
+
|\delta\psi\rangle
+
O(\|\delta z\|^2),
$$

where

$$
|\delta\psi\rangle
=
\sum_\alpha
\delta z_\alpha
|\partial_\alpha\psi\rangle.
$$

This vector \( |\delta\psi\rangle \) is not itself required to be a valid ECG
wavefunction. It is a velocity-like object: it tells us the direction in Hilbert
space in which the ECG wavefunction moves when the parameters are changed.

So we should not think

$$
|\delta\psi\rangle
\in
\mathcal{M}.
$$

Instead, we should think

$$
|\delta\psi\rangle
\in
T_\psi\mathcal{M}.
$$

The tangent space is not the manifold. It is the local linear approximation to
the manifold at the current point.

## Simple Geometric Example

Take a circle in the plane:

$$
\mathcal{M}
=
\{(\cos\theta,\sin\theta)\}.
$$

At angle \(\theta\), a point on the circle is

$$
r(\theta)
=
(\cos\theta,\sin\theta).
$$

If we change the parameter by \(\delta\theta\), the exact new point is

$$
r(\theta+\delta\theta)
=
(\cos(\theta+\delta\theta),\sin(\theta+\delta\theta)).
$$

This point is still exactly on the circle.

The first-order expansion is

$$
r(\theta+\delta\theta)
=
r(\theta)
+
\delta\theta\,r'(\theta)
+
O(\delta\theta^2).
$$

The tangent vector is

$$
\delta\theta\,r'(\theta).
$$

That tangent vector is not itself a point on the circle. It is a vector in the
tangent line. Also, the approximate point

$$
r(\theta)
+
\delta\theta\,r'(\theta)
$$

usually lies slightly off the circle. The error is second order:

$$
O(\delta\theta^2).
$$

The ECG manifold works the same way, except the ambient space is Hilbert space
instead of the plane.

## What TDVP Actually Does

TDVP chooses the parameter velocity

$$
\dot z
$$

so that

$$
|\dot\psi\rangle
=
\sum_\alpha
\dot z_\alpha
|\partial_\alpha\psi\rangle
$$

is the best tangent-space approximation to the exact Schrodinger velocity:

$$
-iH|\psi\rangle.
$$

Then the code updates the parameters:

$$
z(t)
\longrightarrow
z(t+\Delta t).
$$

The actual state after the update is evaluated as

$$
|\psi(z(t+\Delta t))\rangle.
$$

Because this is again built from ECG parameters, it is again an ECG state.

So the guarantee is not that the tangent vector itself is an ECG wavefunction.
The guarantee is:

$$
z(t+\Delta t)
\text{ is a parameter vector}
\quad\Longrightarrow\quad
|\psi(z(t+\Delta t))\rangle
\in
\mathcal{M}.
$$

## One More Practical Condition

For ECGs, the parameters must also remain physically valid. For example, the
Gaussian widths must keep the basis functions normalizable. Roughly speaking,
the real width matrix must remain positive enough:

$$
\mathrm{Re}(A+B)
\succ
0.
$$

If a numerical step makes the Gaussian width invalid, then the parameter vector
no longer describes a proper normalizable ECG basis function. In practice, the
solver must keep time steps small enough, and the implementation must preserve
the required symmetry and positivity conditions on the Gaussian parameters.

## Key Point

The correct statement is:

$$
|\psi(z+\delta z)\rangle
\in
\mathcal{M},
$$

but

$$
\sum_\alpha
\delta z_\alpha
|\partial_\alpha\psi\rangle
\in
T_\psi\mathcal{M},
$$

not necessarily in \(\mathcal{M}\).

The tangent vector is the local direction of motion. The ECG state remains in
the variational manifold because we update and re-evaluate the ECG ansatz
parameters, not because the tangent vector itself is a valid ECG state.

---

## Question 4

Original question:

> In Question 1, why do we project the exact change into the tangent space? The
> tangent space is not a subset of the variational manifold. If we do that,
> don't we need another operator to project the new state back into the
> variational manifold?

Refined question:

> Since \(T_\psi\mathcal{M}\) is not itself the variational manifold
> \(\mathcal{M}\), why is it legitimate for TDVP to project the Schrodinger
> motion onto \(T_\psi\mathcal{M}\)? How does the state remain inside
> \(\mathcal{M}\)?

## Answer

This is exactly the right subtlety.

The important correction is:

> TDVP does not project the finite new state into the tangent space. TDVP
> projects the instantaneous velocity into the tangent space.

The exact Schrodinger equation gives an instantaneous velocity:

$$
|\dot\psi\rangle_{\mathrm{exact}}
=
-iH|\psi\rangle.
$$

This is not yet a new state. It is a direction of motion at the current state.

TDVP replaces this exact velocity by the best velocity that the variational
ansatz can produce:

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
P_T|\dot\psi\rangle_{\mathrm{exact}}
=
-iP_T H|\psi\rangle.
$$

Here

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
\in
T_\psi\mathcal{M}.
$$

That is allowed because velocities of curves on \(\mathcal{M}\) live in
\(T_\psi\mathcal{M}\).

## The State Is Kept on the Manifold by Updating Parameters

The TDVP velocity is written as

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
\sum_\alpha
\dot z_\alpha
|\partial_\alpha\psi\rangle.
$$

This gives an ordinary differential equation for the parameters:

$$
\dot z
=
f(z).
$$

Then we evolve the parameters:

$$
z(t)
\longrightarrow
z(t+\Delta t).
$$

The new variational state is not taken to be

$$
|\psi(t)\rangle
+
\Delta t\,|\dot\psi\rangle_{\mathrm{TDVP}}.
$$

That linearized object can indeed lie off the manifold.

Instead, the new state is evaluated from the ansatz:

$$
|\psi(t+\Delta t)\rangle
=
|\psi(z(t+\Delta t))\rangle.
$$

Because \(z(t+\Delta t)\) is again a set of ECG parameters, the state is again
inside the ECG manifold:

$$
|\psi(z(t+\Delta t))\rangle
\in
\mathcal{M}.
$$

So no extra Hilbert-space projection operator is needed in the ideal TDVP
formulation. The map

$$
z
\mapsto
|\psi(z)\rangle
$$

itself puts the state back on the manifold.

## A Useful Way to Think About It

There are two different objects:

1. The tangent-space velocity:

$$
|\dot\psi\rangle
\in
T_\psi\mathcal{M}.
$$

2. The actual evolved variational state:

$$
|\psi(z(t+\Delta t))\rangle
\in
\mathcal{M}.
$$

The first object is used to decide how the parameters should move. The second
object is the actual state after the parameters are updated.

## Circle Example Again

For a circle,

$$
r(\theta)
=
(\cos\theta,\sin\theta).
$$

The tangent velocity is

$$
\dot r
=
\dot\theta\,r'(\theta).
$$

This vector lies in the tangent line, not on the circle. But if we update the
parameter,

$$
\theta
\longrightarrow
\theta+\Delta\theta,
$$

then the new point is

$$
r(\theta+\Delta\theta)
=
(\cos(\theta+\Delta\theta),\sin(\theta+\Delta\theta)),
$$

which is exactly on the circle.

We do not need to project the tangent vector itself onto the circle. We use the
tangent vector to update the coordinate \(\theta\), and then evaluate the circle
parametrization again.

The ECG case is the same:

$$
z
\longrightarrow
z+\Delta z,
\qquad
|\psi\rangle
\longrightarrow
|\psi(z+\Delta z)\rangle.
$$

## Where Projection Really Happens

The projector \(P_T\) appears only in the velocity equation:

$$
|\dot\psi\rangle
=
-iP_T H|\psi\rangle.
$$

It answers this question:

> Which instantaneous direction inside the tangent space best approximates the
> exact Schrodinger direction?

It does not answer this question:

> Given an arbitrary Hilbert-space vector, how do I project it back onto the
> nonlinear ECG manifold?

That second question would be a different nonlinear fitting problem:

$$
\min_{z'}
\left\|
|\psi(z')\rangle
-
|\varphi\rangle
\right\|^2.
$$

TDVP avoids solving this nonlinear projection at every step. Instead, it
evolves the parameters directly.

## What About Numerical Integration Error?

With a finite time step, an integrator such as Euler or RK4 only approximates
the exact parameter flow:

$$
\dot z
=
f(z).
$$

But after every step, the code still stores ECG parameters and evaluates an ECG
state. So the state remains on the parameterized ECG manifold, provided the
parameters remain physically valid.

The finite-step error is not that the state leaves \(\mathcal{M}\). The
finite-step error is that the computed parameter curve may approximate the true
TDVP parameter curve imperfectly.

## Key Point

The projection is applied to the instantaneous velocity:

$$
-iH|\psi\rangle
\longrightarrow
-iP_T H|\psi\rangle
\in
T_\psi\mathcal{M}.
$$

The actual state remains in the variational manifold because TDVP evolves the
parameters:

$$
z(t)
\longrightarrow
z(t+\Delta t),
$$

and then evaluates

$$
|\psi(z(t+\Delta t))\rangle.
$$

So there is no contradiction: the tangent space is not a subset of the
manifold, but it is the correct space for velocities of curves that lie on the
manifold.

---

## Question 5

Original question:

> I guess the logic is: you minimize the residual between exact velocity and
> first-order velocity. The \(z\) which makes the residual smallest should be
> used to evolve the state in the variational manifold.

Refined question:

> Is TDVP choosing the variational parameter motion by minimizing the
> difference between the exact Schrodinger velocity and the first-order
> variational velocity?

## Answer

Yes, that is the right logic. The only correction is:

> At the current parameter value \(z(t)\), TDVP chooses the best velocity
> \(\dot z(t)\), not the next parameter \(z(t+\Delta t)\) directly.

Then a time integrator uses \(\dot z(t)\) to evolve the parameters.

## Step 1: Current State

At time \(t\), the current variational state is

$$
|\psi(t)\rangle
=
|\psi(z(t))\rangle.
$$

The exact Schrodinger velocity at this same state would be

$$
|\dot\psi\rangle_{\mathrm{exact}}
=
-iH|\psi(z(t))\rangle.
$$

This is the direction the exact wavefunction wants to move.

## Step 2: Variational First-Order Velocity

If the parameters move with velocity \(\dot z\), then to first order

$$
|\psi(z(t+\Delta t))\rangle
=
|\psi(z(t))\rangle
+
\Delta t
\sum_\alpha
\dot z_\alpha
|\partial_\alpha\psi\rangle
+
O(\Delta t^2).
$$

Therefore the variational velocity is

$$
|\dot\psi\rangle_{\mathrm{var}}
=
\sum_\alpha
\dot z_\alpha
|\partial_\alpha\psi\rangle.
$$

Using the tangent map \(T\), this is

$$
|\dot\psi\rangle_{\mathrm{var}}
=
T\dot z.
$$

## Step 3: Minimize the Velocity Error

TDVP chooses \(\dot z\) so that the variational velocity is as close as possible
to the exact velocity:

$$
\min_{\dot z}
\left\|
T\dot z
-
(-iH|\psi\rangle)
\right\|^2.
$$

Equivalently,

$$
\min_{\dot z}
\left\|
T\dot z
+
iH|\psi\rangle
\right\|^2.
$$

Multiplying the residual by \(i\) does not change its norm, so this is the same
as the usual Dirac-Frenkel form:

$$
\min_{\dot z}
\left\|
iT\dot z
-
H|\psi\rangle
\right\|^2.
$$

Since

$$
|\dot\psi\rangle_{\mathrm{var}}
=
T\dot z,
$$

this is also

$$
\min_{\dot z}
\left\|
(i\partial_t-H)|\psi\rangle
\right\|^2.
$$

## Step 4: Solve for \(\dot z\)

The minimizer satisfies the normal equations

$$
iC\dot z
=
g,
$$

where

$$
C_{\alpha\beta}
=
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle,
\qquad
g_\alpha
=
\langle\partial_\alpha\psi|H|\psi\rangle.
$$

If \(C\) is invertible,

$$
\dot z
=
-iC^{-1}g.
$$

So TDVP has produced an ordinary differential equation for the variational
parameters:

$$
\frac{dz}{dt}
=
f(z).
$$

## Step 5: Evolve the Parameters, Then Rebuild the State

For a small time step, a simple Euler update would be

$$
z(t+\Delta t)
\approx
z(t)
+
\Delta t\,\dot z(t).
$$

A higher-order method such as RK4 does this more accurately, but the idea is
the same.

After the parameter update, the new state is evaluated as

$$
|\psi(t+\Delta t)\rangle
=
|\psi(z(t+\Delta t))\rangle.
$$

That state is inside the variational manifold because it is built from ECG
parameters.

## Corrected One-Sentence Logic

Your sentence can be made precise as:

> TDVP minimizes the residual between the exact Schrodinger velocity and the
> first-order velocity available inside the variational manifold. The resulting
> best parameter velocity \(\dot z\) is then used to evolve the parameters
> \(z(t)\), and the updated parameters define the next variational state.

---

## Question 6

Original question:

> I understand the process, but I do not understand why it is right. It feels
> crazy.

Refined question:

> Why is it reasonable to replace the exact Schrodinger velocity by its best
> tangent-space approximation? Is TDVP actually correct, or only an
> approximation?

## Answer

Your reaction is reasonable. TDVP looks strange because it replaces the true
Hilbert-space motion by motion constrained to a much smaller nonlinear
manifold.

The key point is:

> TDVP is not guaranteed to reproduce the exact Schrodinger trajectory. It is
> the locally best possible trajectory inside the chosen variational manifold.

So TDVP is "right" in a conditional sense:

> If we insist that the wavefunction must always remain inside
> \(\mathcal{M}\), then TDVP gives the best instantaneous velocity available at
> each point of \(\mathcal{M}\).

It is not magic. It is a constrained approximation.

## Exact Dynamics Versus Constrained Dynamics

The exact equation wants

$$
|\dot\psi\rangle
=
-iH|\psi\rangle.
$$

But the variational ansatz only allows velocities of the form

$$
|\dot\psi\rangle_{\mathrm{var}}
=
\sum_\alpha
\dot z_\alpha
|\partial_\alpha\psi\rangle
\in
T_\psi\mathcal{M}.
$$

If

$$
-iH|\psi\rangle
\in
T_\psi\mathcal{M},
$$

then TDVP is exact at that point.

But usually

$$
-iH|\psi\rangle
\notin
T_\psi\mathcal{M}.
$$

Then no possible parameter velocity can reproduce the exact motion. Once that
happens, we must choose an approximation. TDVP chooses the closest one:

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
\operatorname*{argmin}_{v\in T_\psi\mathcal{M}}
\|v-(-iH|\psi\rangle)\|^2.
$$

That is why the projection appears:

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
P_T(-iH|\psi\rangle)
=
-iP_T H|\psi\rangle.
$$

So the reason is not that the tangent-space velocity is exactly correct. The
reason is that it is the least-wrong velocity available inside the ansatz.

## Why Local Velocity Matching Makes Sense

For a very small time step \(\Delta t\), the exact state is

$$
|\psi_{\mathrm{exact}}(t+\Delta t)\rangle
=
|\psi(t)\rangle
+
\Delta t\,|\dot\psi\rangle_{\mathrm{exact}}
+
O(\Delta t^2).
$$

A variational trajectory has

$$
|\psi_{\mathrm{var}}(t+\Delta t)\rangle
=
|\psi(t)\rangle
+
\Delta t\,|\dot\psi\rangle_{\mathrm{var}}
+
O(\Delta t^2).
$$

Subtract them:

$$
|\psi_{\mathrm{var}}(t+\Delta t)\rangle
-
|\psi_{\mathrm{exact}}(t+\Delta t)\rangle
=
\Delta t
\left(
|\dot\psi\rangle_{\mathrm{var}}
-
|\dot\psi\rangle_{\mathrm{exact}}
\right)
+
O(\Delta t^2).
$$

To make the one-step error as small as possible to leading order, we should
minimize

$$
\left\|
|\dot\psi\rangle_{\mathrm{var}}
-
|\dot\psi\rangle_{\mathrm{exact}}
\right\|^2.
$$

That is exactly TDVP.

So TDVP is justified because it minimizes the leading-order local error in the
wavefunction.

## Analogy: Constrained Motion on a Surface

Imagine a particle forced to stay on a surface. Suppose the unconstrained
velocity points partly off the surface.

The particle cannot use the off-surface component. The best allowed velocity is
the tangent component.

That does not mean the tangent plane is the surface. It means the tangent plane
contains the allowed instantaneous velocities.

TDVP is the same idea:

$$
\text{exact velocity}
\quad
\longrightarrow
\quad
\text{best allowed tangent velocity}.
$$

Then the parameters are updated, and the state moves along the curved manifold.

## When TDVP Is Exactly Right

TDVP becomes exact if the manifold is invariant under the Hamiltonian flow.
That means:

$$
-iH|\psi\rangle
\in
T_\psi\mathcal{M}
\qquad
\text{for every } |\psi\rangle\in\mathcal{M}.
$$

In that case, the exact Schrodinger velocity is already available inside the
ansatz, so projection changes nothing:

$$
P_T H|\psi\rangle
=
H|\psi\rangle.
$$

Then

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
|\dot\psi\rangle_{\mathrm{exact}}.
$$

This is why Gaussian ansatz methods work extremely well for some quadratic
Hamiltonians: Gaussian states can be closed under the exact dynamics.

## When TDVP Is Only Approximate

If the exact dynamics leaves the manifold, then TDVP is approximate.

The error is the part of the exact velocity that the tangent space cannot
represent:

$$
|\dot\psi\rangle_{\mathrm{exact}}
-
|\dot\psi\rangle_{\mathrm{TDVP}}
=
-i(1-P_T)H|\psi\rangle.
$$

If this residual is small, TDVP is good.

If this residual is large, TDVP can be bad.

That is exactly what the Step-4 report is diagnosing. The ECG basis was good
enough for the ground state, but the real-time trajectory needed tangent
directions that were missing or numerically truncated.

## Why the Method Is Still Useful

TDVP is useful because it turns an impossible problem,

$$
\text{evolve in the full Hilbert space},
$$

into a smaller problem,

$$
\text{evolve the variational parameters } z(t).
$$

Among all possible parameter velocities, it chooses the one with the smallest
instantaneous Schrodinger residual.

This is similar in spirit to the Rayleigh-Ritz variational method for
eigenstates. Rayleigh-Ritz does not give the exact ground state unless the
trial space is complete, but inside the chosen trial space it gives the best
energy approximation.

TDVP is the time-dependent version:

> Inside the chosen tangent space, choose the best dynamical approximation.

## The Important Warning

TDVP is locally optimal, not globally guaranteed.

It minimizes the instantaneous error:

$$
\left\|
|\dot\psi\rangle_{\mathrm{var}}
-
|\dot\psi\rangle_{\mathrm{exact}}
\right\|.
$$

It does not guarantee that the accumulated long-time trajectory remains close
to the exact trajectory.

This is why the verification report found a failure:

$$
\text{energy conserved}
\quad
\text{but}
\quad
\text{norm and recurrence failed}.
$$

The method was still following a projected variational dynamics, but the
projected dynamics was not close enough to the true Schrodinger dynamics after
SVD truncation.

## Short Answer

TDVP is reasonable because, for an infinitesimal time step, matching the exact
velocity as closely as possible minimizes the leading-order error.

It is not guaranteed to be globally correct. It is correct only if the manifold
can represent the exact velocity along the trajectory. Otherwise, it is the
best local approximation available inside the chosen variational ansatz.

---

## Question 7

Original question:

> In Question 2, I guess the matrix language is just a tool. The elements of
> the matrix do not have the physical picture that quantum mechanics describes
> in textbooks.

Refined question:

> Are the matrices \(C\) and \(g\) in the TDVP derivation physical observables,
> or are they just coordinate tools for expressing the tangent-space projection?

## Answer

Yes, the matrix language is mainly a tool.

The real physical/geometric statement is not the matrix equation itself. The
real statement is:

$$
\text{choose the tangent-space velocity closest to}
\quad
-iH|\psi\rangle.
$$

The matrix equation

$$
iC\dot z
=
g
$$

is just the coordinate expression of that projection after we choose a
particular tangent basis

$$
\{|\partial_\alpha\psi\rangle\}.
$$

## What Is Physical and What Is Coordinate-Dependent?

The physical Hilbert-space objects are:

$$
|\psi\rangle,
\qquad
H|\psi\rangle,
\qquad
T_\psi\mathcal{M},
\qquad
P_T H|\psi\rangle,
\qquad
|R\rangle.
$$

These have direct geometric meaning.

But the individual matrix elements

$$
C_{\alpha\beta},
\qquad
g_\alpha,
$$

depend on the chosen coordinates \(z_\alpha\). If we reparameterize the same
manifold using different variables, the numbers in \(C\) and \(g\) change.

So an individual element like

$$
C_{23}
$$

is not usually a physical observable. It is a coordinate-dependent number.

## What \(C\) Means

The matrix

$$
C_{\alpha\beta}
=
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle
$$

is the Gram matrix of tangent vectors.

It tells us how much two parameter directions overlap as wavefunction changes.

The diagonal element

$$
C_{\alpha\alpha}
=
\langle\partial_\alpha\psi|\partial_\alpha\psi\rangle
$$

measures how strongly the wavefunction changes when parameter \(z_\alpha\) is
changed.

The off-diagonal element

$$
C_{\alpha\beta}
=
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle
$$

measures whether two parameter changes produce similar wavefunction changes.

So \(C\) has a geometric meaning:

> \(C\) is the metric tensor on the variational manifold, written in the
> coordinates \(z_\alpha\).

It is not a textbook observable like position, momentum, or energy.

## What \(g\) Means

The vector

$$
g_\alpha
=
\langle\partial_\alpha\psi|H|\psi\rangle
$$

measures how much the Hamiltonian-driven direction

$$
H|\psi\rangle
$$

points along the tangent direction

$$
|\partial_\alpha\psi\rangle.
$$

So \(g_\alpha\) is a kind of projected Hamiltonian drive or generalized force.

But again, a single component \(g_\alpha\) depends on the chosen coordinates.
It is not by itself a physical observable.

## Why the Matrix Equation Is Needed

Because the tangent vectors are generally not orthonormal.

Usually,

$$
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle
\neq
\delta_{\alpha\beta}.
$$

So if we know the overlaps \(g_\alpha\), we cannot simply say

$$
\dot z_\alpha
=
-ig_\alpha.
$$

The tangent directions overlap and mix with each other. The matrix \(C\)
corrects for that nonorthogonality:

$$
\dot z
=
-iC^{-1}g.
$$

This is similar to solving a least-squares problem with a nonorthogonal basis.

## Coordinate-Free Version

The clean coordinate-free formula is

$$
|\dot\psi\rangle
=
-iP_T H|\psi\rangle.
$$

This is the physical/geometric statement.

The coordinate formula is

$$
|\dot\psi\rangle
=
\sum_\alpha
\dot z_\alpha
|\partial_\alpha\psi\rangle,
\qquad
iC\dot z
=
g.
$$

This is the computational statement.

Both describe the same TDVP motion, but the first one is conceptually cleaner.

## Relation to Textbook Quantum Mechanics

In textbook quantum mechanics, matrix elements such as

$$
\langle n|H|m\rangle
$$

often have a direct interpretation as Hamiltonian couplings between basis
states.

Here the situation is different. The vectors

$$
|\partial_\alpha\psi\rangle
$$

are tangent vectors generated by changing variational parameters. They are not
usually normalized physical states with a direct experimental meaning.

Therefore

$$
C_{\alpha\beta}
\quad
\text{and}
\quad
g_\alpha
$$

should mostly be understood as geometric/computational quantities.

The physically meaningful result is the final variational velocity:

$$
|\dot\psi\rangle_{\mathrm{TDVP}}
=
-iP_T H|\psi\rangle,
$$

and the residual:

$$
|R\rangle
=
i|\dot\psi\rangle_{\mathrm{TDVP}}
-
H|\psi\rangle.
$$

## Short Answer

The matrices are not the main physical picture. They are a coordinate language
for doing the tangent-space projection.

The physical picture is:

$$
\text{exact velocity}
\quad
\longrightarrow
\quad
\text{best tangent-space velocity}.
$$

The matrix \(C\) tells us the geometry of the tangent coordinates, and \(g\)
tells us how the Hamiltonian points into those tangent directions. Individual
entries of \(C\) and \(g\) are coordinate-dependent, so they should not be
interpreted as textbook observables.

---

## Question 8

Original question:

> I have some confusion about the size of the state \(|\psi\rangle\) and the
> size of \(T\).

Refined question:

> In the matrix notation
> 
> $$
> |\dot\psi\rangle
> =
> T\dot z,
> \qquad
> C
> =
> T^\dagger T,
> \qquad
> g
> =
> T^\dagger H|\psi\rangle,
> $$
> 
> what are the dimensions of \(|\psi\rangle\), \(T\), \(C\), \(g\), and
> \(\dot z\)?

## Answer

First, \(|\psi\rangle\) is not really a matrix. It is a vector in Hilbert
space. If we represent Hilbert space with \(N_H\) basis functions or grid
points, then we may write it as a column vector:

$$
|\psi\rangle
\sim
\begin{bmatrix}
\psi_1\\
\psi_2\\
\vdots\\
\psi_{N_H}
\end{bmatrix}.
$$

So in a finite representation,

$$
|\psi\rangle
\in
\mathbb{C}^{N_H}.
$$

For a true continuous quantum problem, \(N_H\) is infinite. In that case the
same formulas are formal Hilbert-space formulas, not ordinary finite matrix
multiplications.

## Size of the Parameter Vector

The variational parameters are

$$
z
=
(z_1,z_2,\dots,z_d).
$$

Here \(d\) is the number of variational parameters. Therefore

$$
z
\in
\mathbb{C}^{d},
\qquad
\dot z
\in
\mathbb{C}^{d}.
$$

So \(\dot z\) is a \(d\times 1\) column vector:

$$
\dot z
\sim
\begin{bmatrix}
\dot z_1\\
\dot z_2\\
\vdots\\
\dot z_d
\end{bmatrix}.
$$

## Size of \(T\)

The tangent map \(T\) has columns equal to tangent vectors:

$$
T
=
\begin{bmatrix}
|\partial_1\psi\rangle
&
|\partial_2\psi\rangle
&
\cdots
&
|\partial_d\psi\rangle
\end{bmatrix}.
$$

Each column is a Hilbert-space vector:

$$
|\partial_\alpha\psi\rangle
\in
\mathbb{C}^{N_H}.
$$

There are \(d\) such columns. Therefore, in a finite representation,

$$
T
\in
\mathbb{C}^{N_H\times d}.
$$

So \(T\) maps parameter velocities into wavefunction velocities:

$$
T:
\mathbb{C}^{d}
\longrightarrow
\mathbb{C}^{N_H}.
$$

That is why

$$
|\dot\psi\rangle
=
T\dot z
$$

has size

$$
(N_H\times d)(d\times 1)
=
N_H\times 1.
$$

So the result is again a Hilbert-space vector, as it should be.

## Size of \(C\)

The Gram matrix is

$$
C
=
T^\dagger T.
$$

Since

$$
T^\dagger
\in
\mathbb{C}^{d\times N_H},
\qquad
T
\in
\mathbb{C}^{N_H\times d},
$$

we get

$$
C
\in
\mathbb{C}^{d\times d}.
$$

In components,

$$
C_{\alpha\beta}
=
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle.
$$

So \(C\) is finite as long as the number of variational parameters \(d\) is
finite.

This is important: even if Hilbert space is infinite-dimensional, \(C\) is only
a \(d\times d\) matrix because it lives in parameter/tangent-coordinate space.

## Size of \(g\)

The vector \(g\) is

$$
g
=
T^\dagger H|\psi\rangle.
$$

The sizes are

$$
H|\psi\rangle
\in
\mathbb{C}^{N_H},
$$

and

$$
T^\dagger
\in
\mathbb{C}^{d\times N_H}.
$$

Therefore

$$
g
\in
\mathbb{C}^{d}.
$$

So \(g\) is a \(d\times 1\) vector:

$$
g_\alpha
=
\langle\partial_\alpha\psi|H|\psi\rangle.
$$

## Size Check of the TDVP Equation

The TDVP equation is

$$
iC\dot z
=
g.
$$

Check the dimensions:

$$
C
\in
\mathbb{C}^{d\times d},
\qquad
\dot z
\in
\mathbb{C}^{d\times 1},
$$

so

$$
C\dot z
\in
\mathbb{C}^{d\times 1}.
$$

And

$$
g
\in
\mathbb{C}^{d\times 1}.
$$

Therefore the equation is dimensionally consistent:

$$
iC\dot z
=
g
\qquad
\text{is a } d\times 1 \text{ equation}.
$$

## Size of the Projector \(P_T\)

The tangent-space projector is

$$
P_T
=
TC^{-1}T^\dagger.
$$

Check its size:

$$
T
\in
\mathbb{C}^{N_H\times d},
\qquad
C^{-1}
\in
\mathbb{C}^{d\times d},
\qquad
T^\dagger
\in
\mathbb{C}^{d\times N_H}.
$$

So

$$
P_T
\in
\mathbb{C}^{N_H\times N_H}.
$$

That makes sense: \(P_T\) acts on Hilbert-space vectors.

For example,

$$
P_T H|\psi\rangle
\in
\mathbb{C}^{N_H}.
$$

But in real ECG code, we usually do not build this huge projector explicitly.
Instead, we solve the smaller \(d\times d\) system

$$
iC\dot z
=
g.
$$

## Summary Table

If Hilbert space is represented with size \(N_H\), and the variational ansatz
has \(d\) parameters:

| object | meaning | size |
|---|---|---|
| \(|\psi\rangle\) | wavefunction | \(N_H\times 1\) |
| \(H\) | Hamiltonian operator | \(N_H\times N_H\) |
| \(z\), \(\dot z\) | parameters and parameter velocity | \(d\times 1\) |
| \(T\) | tangent map | \(N_H\times d\) |
| \(T^\dagger\) | maps Hilbert vectors to tangent coordinates | \(d\times N_H\) |
| \(C=T^\dagger T\) | tangent Gram matrix | \(d\times d\) |
| \(g=T^\dagger H|\psi\rangle\) | projected Hamiltonian drive | \(d\times 1\) |
| \(P_T=TC^{-1}T^\dagger\) | Hilbert-space projector | \(N_H\times N_H\) |

## Important ECG Code Point

In the ECG implementation, \(|\psi\rangle\) is not usually stored as a giant
grid vector when building \(C\) and \(g\). The code computes inner products like

$$
\langle\partial_\alpha\psi|\partial_\beta\psi\rangle
$$

and

$$
\langle\partial_\alpha\psi|H|\psi\rangle
$$

analytically using Gaussian integral formulas.

So \(N_H\) may be infinite in principle, but the actual linear algebra we solve
is finite:

$$
iC\dot z
=
g,
\qquad
C\in\mathbb{C}^{d\times d}.
$$

That is one major reason variational methods are useful.

---

## Question 9

Original question:

> Where do the two TDVP conservation laws
> 
> $$
> \frac{d}{dt}\langle\psi|\psi\rangle = 2\,\mathrm{Im}\,\langle\psi|PH|\psi\rangle,
> \qquad
> \frac{d}{dt}\langle\psi|H|\psi\rangle = 2\,\mathrm{Im}\,\langle\psi|HPH|\psi\rangle
> $$
> 
> come from?

Refined question:

> In Section 9.1, the report claims that "direct differentiation" of the
> projected Schrodinger equation
> 
> $$
> |\dot\psi\rangle = -i\,P\,H|\psi\rangle
> $$
> 
> gives the two TDVP conservation laws above for any Hermitian projector \(P\).
> How are these two formulas derived line by line, and why do they immediately
> imply that energy is conserved structurally while norm conservation is only
> conditional?

## Answer

The starting point is the projected flow

$$
|\dot\psi\rangle
=
-iPH|\psi\rangle,
$$

with the assumptions

$$
H^\dagger=H,
\qquad
P^\dagger=P,
\qquad
P^2=P.
$$

We will compute the time derivatives of

$$
\langle\psi|\psi\rangle
\quad\text{and}\quad
\langle\psi|H|\psi\rangle.
$$

## A Single Master Identity

For any time-independent Hermitian operator \(A\),

$$
\frac{d}{dt}\langle\psi|A|\psi\rangle
=
\langle\dot\psi|A|\psi\rangle
+
\langle\psi|A|\dot\psi\rangle.
$$

The two terms are complex conjugates of each other:

$$
\langle\psi|A|\dot\psi\rangle
=
\bigl(\langle\dot\psi|A^\dagger|\psi\rangle\bigr)^{*}
=
\bigl(\langle\dot\psi|A|\psi\rangle\bigr)^{*},
$$

so adding them just doubles the real part:

$$
\frac{d}{dt}\langle\psi|A|\psi\rangle
=
2\,\mathrm{Re}\,\langle\dot\psi|A|\psi\rangle.
\tag{$\ast$}
$$

We also need the bra form of the TDVP velocity. Taking the Hermitian conjugate
of \(|\dot\psi\rangle=-iPH|\psi\rangle\),

$$
\langle\dot\psi|
=
+i\,\langle\psi|H^\dagger P^\dagger
=
+i\,\langle\psi|HP,
$$

using \(H^\dagger=H\) and \(P^\dagger=P\).

A second identity we will reuse: for any complex \(z=a+ib\),

$$
\mathrm{Re}(iz)
=
\mathrm{Re}(ia-b)
=
-b
=
-\mathrm{Im}(z).
\tag{$\dagger$}
$$

## Step 1: The Norm Law

Apply \((\ast)\) with \(A=\mathbb{1}\):

$$
\frac{d}{dt}\langle\psi|\psi\rangle
=
2\,\mathrm{Re}\,\langle\dot\psi|\psi\rangle.
$$

Substitute the bra:

$$
\langle\dot\psi|\psi\rangle
=
(i\langle\psi|HP)|\psi\rangle
=
i\,\langle\psi|HP|\psi\rangle.
$$

Use \((\dagger)\):

$$
\frac{d}{dt}\langle\psi|\psi\rangle
=
-2\,\mathrm{Im}\,\langle\psi|HP|\psi\rangle.
$$

To match the form quoted in the report, take a complex conjugate. For Hermitian
\(H\) and \(P\),

$$
\langle\psi|HP|\psi\rangle^{*}
=
\langle\psi|P^\dagger H^\dagger|\psi\rangle
=
\langle\psi|PH|\psi\rangle,
$$

so

$$
\mathrm{Im}\,\langle\psi|HP|\psi\rangle
=
-\mathrm{Im}\,\langle\psi|PH|\psi\rangle,
$$

and therefore

$$
\boxed{\;
\frac{d}{dt}\langle\psi|\psi\rangle
=
2\,\mathrm{Im}\,\langle\psi|PH|\psi\rangle.
\;}
$$

## Step 2: The Energy Law

Apply \((\ast)\) with \(A=H\):

$$
\frac{d}{dt}\langle\psi|H|\psi\rangle
=
2\,\mathrm{Re}\,\langle\dot\psi|H|\psi\rangle.
$$

Substitute the bra:

$$
\langle\dot\psi|H|\psi\rangle
=
(i\langle\psi|HP)H|\psi\rangle
=
i\,\langle\psi|HPH|\psi\rangle.
$$

Use \((\dagger)\) again:

$$
\frac{d}{dt}\langle\psi|H|\psi\rangle
=
-2\,\mathrm{Im}\,\langle\psi|HPH|\psi\rangle.
$$

The operator \(HPH\) is Hermitian whenever \(H\) and \(P\) are (see Step 3),
so the bilinear form is real and \(\mathrm{Im}\) of it is zero. The sign
therefore does not matter, and the report writes the equivalent form

$$
\boxed{\;
\frac{d}{dt}\langle\psi|H|\psi\rangle
=
2\,\mathrm{Im}\,\langle\psi|HPH|\psi\rangle.
\;}
$$

## Step 3: Energy Is Conserved Structurally

Look at the operator sandwiched in the energy law:

$$
(HPH)^\dagger
=
H^\dagger P^\dagger H^\dagger
=
HPH.
$$

So \(HPH\) is Hermitian, and therefore

$$
\langle\psi|HPH|\psi\rangle
\in
\mathbb{R},
\qquad
\mathrm{Im}\,\langle\psi|HPH|\psi\rangle
=
0.
$$

Hence

$$
\frac{d}{dt}\langle\psi|H|\psi\rangle
=
0
$$

identically, no matter which subspace \(P\) projects onto. Energy conservation
is purely a consequence of \(P\) being a Hermitian projector. It does **not**
require \(P\) to project onto \(T_\psi\mathcal{M}\), nor even for
\(|\psi\rangle\) to lie in \(\mathrm{range}(P)\).

That is why HTML Section 9.4 reports
\(|E(t)-E(0)|/|E(0)|\sim 10^{-13}\) at the round-off floor: the SVD-regularised
projector \(P_{T_r}\) is still Hermitian and idempotent, so the same algebraic
cancellation goes through unchanged.

## Step 4: Norm Conservation Is Only Conditional

The norm law contains the bilinear form

$$
\langle\psi|PH|\psi\rangle.
$$

Now \(PH\) is **not** Hermitian in general:

$$
(PH)^\dagger
=
H^\dagger P^\dagger
=
HP,
$$

and \(HP\neq PH\) unless \([P,H]=0\). So the imaginary part can be nonzero, and
the norm can drift.

There is, however, a clean sufficient condition. If

$$
P|\psi\rangle
=
|\psi\rangle,
\qquad
\text{i.e.}\qquad
|\psi\rangle\in\mathrm{range}(P),
$$

then \(\langle\psi|P=\langle\psi|\) as well, hence

$$
\langle\psi|PH|\psi\rangle
=
\langle\psi|H|\psi\rangle
\in
\mathbb{R},
$$

so \(\mathrm{Im}\,\langle\psi|PH|\psi\rangle=0\) and the norm is preserved.

In short:

> Energy conservation needs only \(P^\dagger=P\). Norm conservation needs the
> additional condition \(|\psi\rangle\in\mathrm{range}(P)\).

## Step 5: Why Exact TDVP Satisfies the Extra Condition

For the full TDVP projector \(P_T\) onto \(T_\psi\mathcal{M}\), the condition
\(|\psi\rangle\in\mathrm{range}(P_T)\) holds automatically in ECG.

The reason is that the linear coefficients \(u_k\) are themselves variational
parameters. Therefore one set of tangent directions is

$$
\frac{\partial|\psi\rangle}{\partial u_k}
=
|\phi_k\rangle,
$$

and the state itself is a linear combination of these tangent vectors:

$$
|\psi\rangle
=
\sum_{k=1}^{K}
u_k|\phi_k\rangle
=
\sum_{k=1}^{K}
u_k\,\partial_{u_k}|\psi\rangle
\in
T_\psi\mathcal{M}.
$$

Hence

$$
|\psi\rangle\in T_\psi\mathcal{M}
=
\mathrm{range}(P_T),
\qquad
P_T|\psi\rangle=|\psi\rangle.
$$

Substituting into Step 4 gives

$$
\frac{d}{dt}\langle\psi|\psi\rangle
=
0
$$

for the exact, un-regularised TDVP flow. Both energy and norm are then
conserved.

## Step 6: Why SVD Truncation Breaks the Norm

The regularised projector \(P_{T_r}\) drops the small singular directions of
the Gram matrix \(C\). Its range

$$
T_r
=
\mathrm{range}(P_{T_r})
\;\subsetneq\;
T_\psi\mathcal{M}
$$

is strictly smaller. After truncation, \(|\psi\rangle\) still lies in
\(T_\psi\mathcal{M}\) but generically **not** in \(T_r\). The condition
\(P_{T_r}|\psi\rangle=|\psi\rangle\) fails:

$$
|\psi\rangle
=
\underbrace{P_{T_r}|\psi\rangle}_{=\,|\psi_T\rangle}
+
\underbrace{(\mathbb{1}-P_{T_r})|\psi\rangle}_{=\,|\psi_\perp\rangle\,\neq\,0}.
$$

Plugging into the norm law,

$$
\frac{d}{dt}\langle\psi|\psi\rangle
=
2\,\mathrm{Im}\,\langle\psi|P_{T_r} H|\psi\rangle
=
2\,\mathrm{Im}\,\langle\psi_T|H|\psi\rangle
=
2\,\mathrm{Im}\,\langle\psi_T|H|\psi_\perp\rangle,
$$

where the last equality uses
\(\langle\psi_T|H|\psi_T\rangle\in\mathbb{R}\). That single off-diagonal matrix
element generically has a nonzero imaginary part and produces the monotonic
norm leak quantified in HTML Section 9.3.

## Compact Summary

Both conservation laws come from the same one-line calculation. For any
time-independent Hermitian \(A\),

$$
\frac{d}{dt}\langle\psi|A|\psi\rangle
=
2\,\mathrm{Re}\,\langle\dot\psi|A|\psi\rangle
=
-2\,\mathrm{Im}\,\langle\psi|HPA|\psi\rangle
=
2\,\mathrm{Im}\,\langle\psi|APH|\psi\rangle.
$$

Specialising:

| \(A\) | conservation law | requires |
|---|---|---|
| \(\mathbb{1}\) | \(\frac{d}{dt}\langle\psi|\psi\rangle=2\,\mathrm{Im}\,\langle\psi|PH|\psi\rangle\) | \(|\psi\rangle\in\mathrm{range}(P)\) for it to vanish |
| \(H\) | \(\frac{d}{dt}\langle H\rangle=2\,\mathrm{Im}\,\langle\psi|HPH|\psi\rangle\) | only \(P^\dagger=P\); vanishes identically |

This is the precise reason the Step-4 verification sees energy at the round-off
floor while the norm leaks to \(0.47\) over \(T=2\pi\): SVD truncation keeps
\(P_{T_r}\) Hermitian (so energy stays exact) but shrinks
\(\mathrm{range}(P_{T_r})\) so that \(|\psi\rangle\) escapes it (so norm
drifts).

---

## Question 10

Original question:

> Now I have an overall understanding of real-time TDVP. So I wonder whether
> real-time TDVP can describe the dynamics of the system. In other words, which
> quantities should we use to describe the dynamics? The wavefunction?

Refined question:

> Once real-time TDVP has produced a parameter trajectory \(z(t)\) and hence a
> variational state \(|\psi(z(t))\rangle\), what quantities actually represent
> the "dynamics of the system"? Is it the wavefunction itself, or some derived
> object?

## Answer

Short answer: yes, real-time TDVP describes the dynamics &mdash; but the
wavefunction itself is rarely the quantity you actually report. You report
**observables** computed from it.

## What "the Dynamics" Means in Practice

The wavefunction \(|\psi(t)\rangle\) contains everything, but it is a vector in
an infinite-dimensional Hilbert space. You cannot plot it, measure it, or
compare it directly with experiment. So in practice, "describing the dynamics"
means tracking time-dependent expectation values:

$$
\langle A\rangle(t)
=
\langle\psi(t)|A|\psi(t)\rangle,
$$

for a chosen set of operators \(A\). The choice of \(A\) depends on what
physics you are after.

## What TDVP Actually Delivers

TDVP gives you two equivalent objects:

1. The parameter trajectory

$$
z(t)
=
(u(t),A(t),B(t),R(t)),
$$

which is what the code stores at every step.

2. The wavefunction

$$
|\psi(z(t))\rangle,
$$

reconstructed on demand from the parameters.

Once you have either, every observable falls out:

$$
\langle A\rangle(t)
=
\langle\psi(z(t))|A|\psi(z(t))\rangle.
$$

So the wavefunction is a **container**, not the answer.

## Which Observables Matter for MBDL

For the Science 389 paper this project targets, the three experimental
signatures of many-body dynamical localization are explicit observables:

- **Momentum distribution** \(n(k,t)=\langle\psi(t)|\,\hat n_k\,|\psi(t)\rangle\)
  &mdash; freezes after \(\sim 10\) kicks if MBDL holds.
- **Energy and entropy** \(\langle H\rangle(t)\), \(S(t)\) &mdash; saturate
  instead of growing linearly.
- **First-order correlation**
  \(G^{(1)}(z,t)=\langle\psi(t)|\,\hat\psi^\dagger(z)\hat\psi(0)\,|\psi(t)\rangle\)
  &mdash; has a distinct decay shape.

These are the things to plot. None of them is the wavefunction.

## Diagnostic Quantities (Not Physics, But Essential)

Separate from observables, you also track quantities that test whether TDVP is
trustworthy at each time:

- \(\langle\psi|\psi\rangle(t)\) &mdash; should stay 1 (Step 4 saw it leak to
  0.47).
- \(\langle H\rangle(t)\) &mdash; conserved exactly under any Hermitian projector
  (Step 4 confirmed at \(10^{-13}\)).
- TDVP residual \(\|(1-P_T)H|\psi\rangle\|\) &mdash; if this is large, the
  manifold is too small.
- Gram-matrix conditioning \(\sigma_1/\sigma_d\) &mdash; a large condition
  number means the SVD cutoff is throwing away physical tangent directions.

The split is important:

> Observables answer the **physics** question. Diagnostics answer the
> **"should I believe this run?"** question.

The Step-4 report is exactly the second kind.

## So Can TDVP Describe the Dynamics?

Yes, but with the caveat already discussed in Question 6: TDVP is locally
optimal, not globally guaranteed. It gives the best dynamics inside the chosen
variational manifold. If the true trajectory leaves that manifold &mdash; or if
numerical truncation removes the directions the trajectory needs, as in Step 4
&mdash; the observables you compute from \(|\psi(z(t))\rangle\) will deviate
from the exact ones, even when energy is conserved to machine precision.

## Short Answer

The dynamics are described by:

$$
\text{parameter trajectory }
z(t)
\;\Longrightarrow\;
\text{wavefunction }
|\psi(z(t))\rangle
\;\Longrightarrow\;
\text{observables }
\langle A\rangle(t).
$$

The wavefunction is the bridge, not the destination. What you report &mdash;
and what you compare with experiment &mdash; are the time-dependent expectation
values of physically meaningful operators.

---

## Question 11

Original question:

> In real-time TDVP there is no longer any reason to keep the norm constant
> (Question 9 showed exact TDVP conserves it but SVD-truncated TDVP does not).
> Can we still use real-time TDVP to compute \(\langle x\rangle,\langle
> x^2\rangle,\langle p\rangle,\langle p^2\rangle\)? What is the right
> mathematical procedure when the norm is leaking?

Refined question:

> (a) How are the four moments \(\langle x\rangle,\langle x^2\rangle,\langle
> p\rangle,\langle p^2\rangle\) computed from an ECG state? (b) Under
> SVD-truncated TDVP where \(\frac{d}{dt}\langle\psi|\psi\rangle\neq 0\), how
> do raw \(\langle\psi|A|\psi\rangle\) and normalised
> \(\langle A\rangle = \langle\psi|A|\psi\rangle/\langle\psi|\psi\rangle\)
> evolve? Which one represents the physical dynamics, and what does the
> harmonic trap predict?

## Answer

The dynamics are well defined under TDVP even when the norm is not conserved.
The trick is to track **two parallel quantities** at every time:

$$
\langle A\rangle_{\mathrm{raw}}(t)
=
\langle\psi(t)|A|\psi(t)\rangle,
\qquad
\langle A\rangle_{\mathrm{norm}}(t)
=
\frac{\langle\psi(t)|A|\psi(t)\rangle}{\langle\psi(t)|\psi(t)\rangle}.
$$

The raw form is what the conservation laws of Question 9 act on. The
normalised form is what compares with experiment and with reference solvers.
We will derive how each of \(A\in\{x,p,x^2,p^2\}\) is evaluated for an ECG
state, and how each form evolves under the projected flow
\(|\dot\psi\rangle=-iPH|\psi\rangle\).

## Part A: Evaluating the Moments for an ECG State

Write the variational state as a superposition of complex Gaussians. For the
single-particle case (the verification harness Step 4),

$$
|\psi\rangle
=
\sum_{k=1}^{K}
u_k|\phi_k\rangle,
\qquad
\phi_k(x)
\propto
\exp\bigl(-\alpha_k x^2 + \beta_k x\bigr),
$$

with complex \(\alpha_k=A_k+B_k\) (real part positive for normalisability) and
\(\beta_k=2 R_k B_k\). Pair-overlap quantities of the form
\(M_G^{(ij)}=\int\phi_i^{*}\phi_j\,dx\), the centred mean
\(\mu^{(ij)}\), and the joint width \(\alpha^{(ij)}=\alpha_i^{*}+\alpha_j\)
are all standard Gaussian integrals; they are precomputed in
`PairCache` and exposed by `c.M_G`, `c.mu(0)`, `c.K_Mj(0,0)`, `c.g_Mj(0)` in
the ECG code (`src/realtime_tdvp.cpp` lines 326&ndash;355).

### \(\langle x\rangle_{\mathrm{raw}}\) and \(\langle p\rangle_{\mathrm{raw}}\)

For a single basis pair \((i,j)\),

$$
\langle\phi_i|\hat x|\phi_j\rangle
=
M_G^{(ij)}\,\mu^{(ij)},
\qquad
\langle\phi_i|\hat p|\phi_j\rangle
=
\langle\phi_i|-i\partial_x|\phi_j\rangle
=
-i\,M_G^{(ij)}\bigl(-2\alpha_j\mu^{(ij)}+\beta_j\bigr).
$$

Summing with the variational coefficients,

$$
\langle\psi|\hat x|\psi\rangle_{\mathrm{raw}}
=
\sum_{i,j}u_i^{*}u_j\,M_G^{(ij)}\mu^{(ij)},
$$

$$
\langle\psi|\hat p|\psi\rangle_{\mathrm{raw}}
=
\sum_{i,j}u_i^{*}u_j\,(-i)M_G^{(ij)}\bigl(-2\alpha_j\mu^{(ij)}+\beta_j\bigr).
$$

These are exactly `amp_x` and `amp_p` in `sample_observables`.

### \(\langle x^2\rangle_{\mathrm{raw}}\) and \(\langle p^2\rangle_{\mathrm{raw}}\)

Rather than a separate Gaussian moment integral, recognise that the
Hamiltonian functionals already isolate these moments:

$$
V_{\mathrm{ho}}
=
\langle\psi|\tfrac{1}{2}m\omega^2 \hat x^2|\psi\rangle
\;\Longrightarrow\;
\langle\psi|\hat x^2|\psi\rangle_{\mathrm{raw}}
=
\frac{2\,V_{\mathrm{ho}}}{m\omega^2},
$$

$$
T_{\mathrm{kin}}
=
\langle\psi|\tfrac{\hat p^2}{2m}|\psi\rangle
\;\Longrightarrow\;
\langle\psi|\hat p^2|\psi\rangle_{\mathrm{raw}}
=
2m\,T_{\mathrm{kin}}.
$$

So \(\langle x^2\rangle_{\mathrm{raw}}\) and \(\langle p^2\rangle_{\mathrm{raw}}\)
are byproducts of the same `Harmonic_functional` and `kinetic_energy_functional`
calls used for the energy. No new integrals are required.

### Normalised forms

Divide each by \(S=\langle\psi|\psi\rangle\):

$$
\langle x\rangle
=
\frac{\langle\psi|\hat x|\psi\rangle_{\mathrm{raw}}}{S},
\qquad
\langle x^2\rangle
=
\frac{2 V_{\mathrm{ho}}}{m\omega^2 S},
\qquad
\langle p^2\rangle
=
\frac{2 m T_{\mathrm{kin}}}{S}.
$$

This single divide is the entire difference between the two forms. The cost
of computing both is one additional double precision division per trace step.

## Part B: Time Evolution Under the Projected Flow

Now ask how \(\langle A\rangle_{\mathrm{raw}}\) and \(\langle A\rangle_{\mathrm{norm}}\)
change with time when \(|\dot\psi\rangle=-iPH|\psi\rangle\) for a Hermitian
projector \(P\). Take \(A\) Hermitian and time-independent.

### Raw form

Repeat the master identity from Question 9 with \(A\) general:

$$
\frac{d}{dt}\langle\psi|A|\psi\rangle_{\mathrm{raw}}
=
2\,\mathrm{Re}\,\langle\dot\psi|A|\psi\rangle
=
2\,\mathrm{Re}\bigl(i\langle\psi|HPA|\psi\rangle\bigr)
=
-2\,\mathrm{Im}\,\langle\psi|HPA|\psi\rangle.
\tag{$\star$}
$$

This is the **modified Ehrenfest relation** for projected dynamics.

### Recovery of standard Ehrenfest when P = 1

When the projector is the identity (no truncation; ansatz is full Hilbert
space), \(HPA=HA\) and

$$
-2\,\mathrm{Im}\,\langle\psi|HA|\psi\rangle
=
i\,\langle\psi|[H,A]|\psi\rangle,
$$

using \(\langle\psi|HA|\psi\rangle-\langle\psi|AH|\psi\rangle = \langle\psi|[H,A]|\psi\rangle\)
and the identity \(2i\,\mathrm{Im}(z)=z-z^{*}\). So $(\star)$ collapses to the
textbook

$$
\frac{d}{dt}\langle A\rangle
=
i\,\langle[H,A]\rangle.
$$

For the harmonic Hamiltonian \(H=\hat p^2/2m+\tfrac{1}{2}m\omega^2\hat x^2\),
the standard commutators give

$$
\frac{d\langle x\rangle}{dt}
=
\frac{\langle p\rangle}{m},
\qquad
\frac{d\langle p\rangle}{dt}
=
-m\omega^2\langle x\rangle,
$$

so \(\langle x\rangle(t)\) and \(\langle p\rangle(t)\) execute simple harmonic
oscillation at frequency \(\omega\). For the second moments,

$$
\frac{d\langle x^2\rangle}{dt}
=
\frac{\langle\hat x\hat p+\hat p\hat x\rangle}{m},
\qquad
\frac{d\langle p^2\rangle}{dt}
=
-m\omega^2\langle\hat x\hat p+\hat p\hat x\rangle,
$$

which gives the conserved energy combination

$$
\frac{1}{2m}\langle p^2\rangle
+
\frac{1}{2}m\omega^2\langle x^2\rangle
=
\langle H\rangle
=
\text{const}.
$$

These are the predictions any TDVP run must reproduce in the limit of an
ansatz rich enough to contain the exact trajectory.

### What changes when P is a non-trivial projector

For projected TDVP \((\star)\) becomes

$$
\frac{d}{dt}\langle\psi|A|\psi\rangle_{\mathrm{raw}}
=
-2\,\mathrm{Im}\,\langle\psi|HPA|\psi\rangle
\;\neq\;
i\,\langle\psi|[H,A]|\psi\rangle
\quad
\text{in general},
$$

because \(HPA-APH\neq H A-AH\) once \(P\neq\mathbb{1}\). The "exact" Ehrenfest
relation gets a residual correction

$$
\frac{d}{dt}\langle A\rangle_{\mathrm{raw}}
-
i\langle[H,A]\rangle_{\mathrm{raw}}
=
2\,\mathrm{Im}\,\langle\psi|H(\mathbb{1}-P)A|\psi\rangle.
$$

This vanishes when \(P|\psi\rangle=|\psi\rangle\) **and** the right-hand side
operator acts inside the tangent range. For the SVD-truncated projector
\(P_{T_r}\) of HTML &sect;9.2, neither holds exactly and the standard Ehrenfest
identities receive a small residual.

### Special case A = H

This is exactly the energy law of Question 9. \(HPH\) is Hermitian, so

$$
\frac{d}{dt}\langle\psi|H|\psi\rangle_{\mathrm{raw}}
=
-2\,\mathrm{Im}\,\langle\psi|HPH|\psi\rangle
=
0
$$

identically &mdash; energy is conserved structurally regardless of the
projector. This is the basis of the diagnostic mentioned in the
implementation plan: even with the norm leaking, the raw quantity
\(\langle\psi|H|\psi\rangle = \langle H\rangle_{\mathrm{norm}}\cdot S\)
should be flat to round-off.

## Part C: Time Evolution of the Normalised Moments

Apply the quotient rule. With \(S=\langle\psi|\psi\rangle\) and
\(N_A=\langle\psi|A|\psi\rangle\),

$$
\frac{d}{dt}\langle A\rangle_{\mathrm{norm}}
=
\frac{1}{S}\frac{dN_A}{dt}
-
\frac{N_A}{S^2}\frac{dS}{dt}.
$$

Using $(\star)$ with \(A\) and with \(A=\mathbb{1}\) (Question 9 norm law),

$$
\frac{dN_A}{dt}
=
-2\,\mathrm{Im}\,\langle\psi|HPA|\psi\rangle,
\qquad
\frac{dS}{dt}
=
-2\,\mathrm{Im}\,\langle\psi|HP|\psi\rangle.
$$

Therefore

$$
\boxed{\;
\frac{d}{dt}\langle A\rangle_{\mathrm{norm}}
=
\frac{-2}{S}\Bigl[
\mathrm{Im}\,\langle\psi|HPA|\psi\rangle
\;-\;
\langle A\rangle_{\mathrm{norm}}\,\mathrm{Im}\,\langle\psi|HP|\psi\rangle
\Bigr].
\;}
$$

The two imaginary parts are coupled. Two regimes are immediate:

1. **Exact TDVP** (\(P|\psi\rangle=|\psi\rangle\)). Then
   \(\mathrm{Im}\,\langle\psi|HP|\psi\rangle=\mathrm{Im}\,\langle\psi|H|\psi\rangle=0\)
   and the norm-correction term drops, while the first term reduces to
   \(i\langle[H,A]\rangle/S\), giving back the textbook Ehrenfest equation
   for normalised expectation values. Norm and dynamics decouple cleanly.

2. **SVD-truncated TDVP** (\(P_{T_r}|\psi\rangle\neq|\psi\rangle\)). Both
   imaginary parts are generically nonzero. The norm leak (Question 9 step 6)
   feeds back into the moment evolution through the second term. The
   normalised moments are still the physically meaningful ones, but their
   trajectory is no longer the unprojected Ehrenfest trajectory; it is
   shifted by a residual proportional to the imaginary off-diagonal matrix
   element \(\mathrm{Im}\,\langle\psi_T|H|\psi_\perp\rangle\) of HTML &sect;9.3.

## Part D: Concrete Predictions for Step 4 (Harmonic Trap, \(N=1\))

The Step-4 evolution Hamiltonian is \(H=\hat p^2/2m+\tfrac{1}{2}m\omega^2\hat x^2\)
(cos-quench off). With \(m=\omega=1\) and the seed centred at the origin:

- **\(\langle x\rangle_{\mathrm{norm}}(t),\langle p\rangle_{\mathrm{norm}}(t)\)**:
  exact dynamics says these stay at zero by parity. Step 4 already saw a small
  symmetry breaking (HTML &sect;6, &sect;9.5) seeded by the chirp residue plus
  the projector residual; expect the same in the normalised moments.

- **\(\langle x^2\rangle_{\mathrm{norm}}+\langle p^2\rangle_{\mathrm{norm}}/(m\omega^2)
  =2\langle H\rangle_{\mathrm{norm}}/(m\omega^2)\)**: should be flat in time.
  But \(\langle H\rangle_{\mathrm{norm}}\) is **not** structurally conserved
  &mdash; only \(\langle\psi|H|\psi\rangle_{\mathrm{raw}}\) is. Since
  \(\langle H\rangle_{\mathrm{norm}}=\langle\psi|H|\psi\rangle_{\mathrm{raw}}/S\)
  and \(S\) drifts, the normalised energy and hence the conserved combination
  drift along with the norm, by exactly the same factor \(1/S(t)\).

- **Raw moments**: \(\langle x^2\rangle_{\mathrm{raw}}\) and
  \(\langle p^2\rangle_{\mathrm{raw}}\) decay together with \(S(t)\), so their
  ratio holds the physical content. \(\langle\psi|H|\psi\rangle_{\mathrm{raw}}
  =\langle H\rangle_{\mathrm{norm}}\cdot S\) should be conserved to round-off
  &mdash; this is the cleanest diagnostic of whether the projector is
  Hermitian (it is, by construction) and whether the integrator is faithful.

## Part E: Practical Recipe

For each trace sample at time \(t\), compute and store:

1. \(S(t)=\langle\psi|\psi\rangle\) (already done as `norm`).
2. Raw moments:
   $$
   \langle x\rangle_{\mathrm{raw}},\quad
   \langle p\rangle_{\mathrm{raw}},\quad
   \langle x^2\rangle_{\mathrm{raw}}=\frac{2 V_{\mathrm{ho}}}{m\omega^2},\quad
   \langle p^2\rangle_{\mathrm{raw}}=2m T_{\mathrm{kin}}.
   $$
3. Normalised moments: divide each raw value by \(S(t)\).
4. Diagnostic: \(\langle\psi|H|\psi\rangle_{\mathrm{raw}}=E_{\mathrm{norm}}\cdot S\).
   Plot \(|\Delta(E_{\mathrm{norm}}\cdot S)|\) over time &mdash; this is the
   Question 9 conservation law made numerical.

For comparison with grid / DVR reference solvers (which propagate a unitary
flow with \(S=1\)), use only the **normalised** form. Differences
\(\Delta\langle x^2\rangle, \Delta\langle p^2\rangle\) between ECG and
reference quantify how the SVD truncation distorts the physical observables,
separately from the trivial bookkeeping effect of the norm leak.

## Short Answer

Yes &mdash; real-time TDVP can compute all four moments under non-conserved
norm. The recipe:

1. **Compute raw**: \(\langle\psi|A|\psi\rangle\), with \(\langle x^2\rangle_{\mathrm{raw}}\)
   and \(\langle p^2\rangle_{\mathrm{raw}}\) extracted directly from the
   harmonic and kinetic functionals.
2. **Compute normalised**: divide by \(S(t)\). This is what to compare with
   experiment and with reference solvers.
3. **Track both**: their ratio is a direct view of the norm-leak history; the
   raw \(\langle\psi|H|\psi\rangle\) tests the structural conservation law of
   Question 9 numerically.

The dynamics of the normalised moments are governed by a modified Ehrenfest
equation (Part C) that reduces to the textbook form for exact TDVP and picks
up a norm-leak correction for SVD-truncated TDVP.
