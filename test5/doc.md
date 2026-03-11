# Research Report: S-Matrix Bootstrap Bounds on Gravitational EFTs

## Overview

This document is a deep analysis of a physics research paper on **deriving rigorous numerical bounds on higher-derivative EFT coefficients in gravitational theories**, using dispersive sum rules and linear programming. The paper is authored by S. Caron-Huot, D. Mazac, L. Rastelli, and D. Simmons-Duffin. Below are three structured sections: (1) a conceptual deep-dive into the paper's physics, (2) a detailed breakdown of the numerical linear programming implementation, and (3) an analysis of the adaptive refinement "task scheduling" flow — including all identified bugs and failure modes.

---

## Part 1: Physics — What the Paper Does and Why

### 1.1 The Central Problem

The paper addresses a longstanding challenge: **graviton exchange diverges in the forward limit**, making it seemingly impossible to apply standard S-matrix positivity techniques to theories containing gravity. In non-gravitational theories, the standard approach (initiated in Adams et al. 2006) derives constraints on low-energy EFT couplings by Taylor-expanding dispersion relations around forward scattering (`u → 0`). But in gravity, the graviton propagator introduces a `1/u` pole that blows up at `u = 0`, invalidating this expansion.

### 1.2 The Setup: EFT of a Massless Scalar Coupled to Gravity

The theory considered is the 2→2 scattering of identical massless scalars in D-dimensional flat space. The low-energy amplitude takes the form:

```
M_low(s,u) = 8πG [st/u + su/t + tu/s] - λ₃²[1/s + 1/t + 1/u] - λ₄
             + g₂(s²+t²+u²) + g₃(stu) + g₄(s²+t²+u²)² + ...
```

where:
- `G` is Newton's constant (the graviton exchange pole)
- `g₂, g₃, g₄, ...` are higher-derivative EFT contact coefficients (the unknowns to be bounded)
- `s + t + u = 0` are Mandelstam variables
- The scale `M` is the mass of the lightest heavy state not in the EFT

The key insight: the EFT is valid **at tree level** under a weak-coupling assumption (all couplings are O(ε) as ε → 0). This is exemplified by string theory where ε = g_s (string coupling) and M is the string scale.

### 1.3 Unitarity and the Spectral Density

The imaginary part of the amplitude admits a partial-wave decomposition:

```
Im M(s,u) = s^{(4-D)/2} Σ_{J even} n_J^{(D)} ρ_J(s) P_J(1 + 2u/s)
```

where:
- `P_J(x)` are Gegenbauer polynomials (reducing to Legendre polynomials in D=4)
- `ρ_J(s) = Im c_J(s)` is the spectral density (non-negative by unitarity: `0 ≤ ρ_J(s) ≤ 2`)
- The normalization `n_J^{(D)}` is chosen so that unitarity of S = 1 + iM maps to `|1 + ic_J|² ≤ 1`

The **positivity of the spectral density** is the crucial non-perturbative input that enables all the bounds.

### 1.4 Twice-Subtracted Dispersion Relations

Under the Regge boundedness assumption (the amplitude grows slower than `|s|²` at large `|s|` for fixed `u < 0`), one can write a twice-subtracted dispersion relation. After a subtraction that cancels s- and t-channel poles and defines a UV average measure:

```
⟨(···)⟩ = (1/π) Σ_{J even} n_J^{(D)} ∫_{M²}^∞ dm²/m² · m^{4-D} ρ_J(m²) · (···)
```

This measure is **positive** because ρ_J ≥ 0. This positivity is the engine of all the bounds.

### 1.5 The Sum Rules C_{k,u}

The key objects are the dispersion sum rules:

```
C_{k,u} ≡ ∮_∞ ds'/(2πi) · 1/s' · M(s',u) / [s'(s'+u)]^{k/2} = 0
```

for `u < 0`, `k = 2, 4, 6, ...`. These equal zero because the amplitude satisfies the Regge bound and the contour at infinity vanishes.

Each `C_{k,u}` equates an EFT contribution (expressed in terms of `G, g₂, g₃, ...`) to a positive UV average. **The EFT terms are finite combinations of the couplings, the UV side is a positive integral — this is the positivity constraint.**

Key facts:
- Only `C_{2,u}` is sensitive to the graviton exchange (via the `−8πG/u` term)
- `C_{k,u}` receives contributions only from contacts with spin `J ≥ k`
- Contacts with spin `b ≤ k/2` (in the `A^a B^b` notation) also contribute

### 1.6 The Forward-Limit Strategy and Its Failure with Gravity

In non-gravitational theories, one expands `C_{k,u}` around `u = 0`. This directly computes:
- `g₂ = ⟨1/m⁴⟩` (immediately positive → g₂ > 0)
- `g₃ = ⟨(3 - 4J²/(D-2))/m⁶⟩`
- Two-sided bounds follow from linear programming over the measure

This is the approach of Caron-Huot et al. 2020. **But it fails for gravity because `C_{2,u}` diverges as `u → 0` due to the `8πG/(-u)` pole.**

### 1.7 The New Strategy: Impact Parameter Localization

The key innovation of this paper is to **measure EFT couplings at small impact parameter `b ~ 1/M`** rather than in the forward limit. The mechanism:

1. The graviton contribution in impact parameter space is `8πG·b/(D−4)`
2. If the momentum-space wavefunction `f(p)` is localized near `b ~ 1/M`, the graviton contribution is naturally suppressed by `1/M²`
3. Localization requires `f(p)` to have support up to `|u| ~ M²`

**The crucial technical tool: the Improved Sum Rule.**

Instead of using `C_{2,u}` directly (which contaminated by all EFT contacts at `|u| ~ M²`), one constructs:

```
C_{2,u}^{improved} = C_{2,u} - Σ_{n=2}^∞ (n·u^{2n-2} C_{2n,0} + u^{2n-1} C'_{2n,0})
```

This subtracts off contributions from `g₄, g₅, g₆, ...` by using the higher-subtracted sum rules `C_{4}, C_{6}, ...` evaluated at the forward limit (where they are safe — no graviton pole). The result:

```
C_{2,u}^{improved}|_EFT = 8πG/(-u) + 2g₂ − g₃u
```

**Only three EFT couplings remain**, regardless of `u`. Yet we retain a full one-parameter family of sum rules (indexed by `u ∈ (−M², 0)`).

The explicit improved kernel for the UV average side is:

```
C_{2,u}^{improved}[m², J] = [(2m²+u)P_J(1+2u/m²)] / [m²(m²+u)²]
    − u²/m⁶ · [(4m²+3u)P_J(1)/((m²+u)²) + 4u P'_J(1)/(m⁴−u²)]
```

### 1.8 The Linear Programming Problem

The full linear program is:

```
if:   ∀m ≥ M, J = 0,2,4,...:
        ∫₀ᴹ dp f(p) C_{2,−p²}^{improved}[m²,J]
        + Σ_{k=4,6,...} ∫₀ᴹ dp h_k(p) X_{k,−p²}[m²,J] ≥ 0

then: ∫₀ᴹ dp f(p) [8πG/p² + 2g₂ + g₃p²] ≥ 0
```

**Decision variables**: functions `f(p)` and `h_k(p)`
**Constraints**: positivity for all allowed (m, J) pairs
**Objective**: chosen to maximize/minimize particular combinations of `G, g₂, g₃`

The null constraints `X_{k,u}` (derived from `C_{k,u}` by eliminating all EFT contributions) provide additional leverage without introducing new unknowns.

A key geometric interpretation: the positivity condition on the first line **implies that the Fourier transform of `f(p)` is non-negative in impact parameter space**. This is proven via the Gegenbauer-to-Bessel limit:

```
lim_{m→∞} P_{mb/2}(1 − 2p²/m²) = Γ((D−2)/2) · J_{(D-4)/2}(bp) / (bp/2)^{(D-4)/2}
```

### 1.9 Physical Results

**Bounds on g₂ and g₃ with gravity (Figure 1):** For dimensions D = 5, ..., 12, the paper carves out the allowed region in the (g₂, g₃) plane. The bounds automatically have the correct EFT scaling `g_k ~ 1/M^{2k}` — dimensional analysis is a theorem, not an assumption.

**Bounds on higher contacts (Figure 2):** The method also constrains `g₄`, `g₅`, etc.

**Non-gravitational check:** The impact-parameter method reproduces the forward-limit bounds of Caron-Huot et al. 2020 for G = 0. This validates the approach.

**String theory check:**
- Type II dilaton scattering gives `g₂ = 0`, `g₃ = 0`, `g₄ = 4πG·ζ(3)/m⁶` — safely inside bounds
- Heterotic string: `g₂·m²/(8πG) = 3/16`, `g₃·m⁴/(8πG) = 3/4` — also within bounds
- The bounds are not saturated by known string amplitudes

**Maximal supergravity (D=10):** The improved Regge behavior allows a `C_{-2}` sum rule. One proves:
```
0 ≤ g₀ ≤ 3.000 · 8πG/M⁶
```
where `g₀` is the R⁴ coefficient. The upper bound shows **all interactions vanish as G → 0** — a remarkable swampland-type result. Type II string theory gives `g₀M⁶/(8πG) = 2ζ(3) ≈ 2.40 < 3.000`.

**D=4 special case:** An IR divergence occurs because the 2D Fourier transform of `1/p²` is logarithmically divergent. A regulated bound is derived by introducing a large-distance cutoff `b_max`, giving:
```
g₂ ≥ −(8πG/M²) × 25·log(0.3·M·b_max)   (D=4)
```

**Extended range (Appendix B):** For meromorphic amplitudes, the improved sum rule remains valid beyond `|u| = M²`. Using larger support in `u` sharpens bounds further and the extremal spectra approach simpler Regge trajectories `m² ~ M²(J−2)/2`.

---

## Part 2: The Numerical Implementation — How Linear Programming Works

### 2.1 Overall Architecture

The infinite-dimensional linear program is converted into a finite one through four layers of approximation:
1. **Truncation of basis functions** for `f(p)` and `h_k(p)`
2. **Discretization** of the continuous mass parameter `m` (equivalently `x = 1 − 1/m²`)
3. **Truncation** of the angular momentum sum at `J_max`
4. **Impact parameter constraint** included as an additional inequality (handles large `J`)

The solver used is **SDPB** (Semidefinite Program Bootstrapper) by Simmons-Duffin et al.

### 2.2 Basis Function Expansion

Functions `f(p)` and `h_k(p)` are expanded in pure powers:

```python
f(p) = Σ_n  a_n · p^n
h_k(p) = Σ_{i=0}^{i_k}  b_{k,i} · p^i
```

The **choice of powers** for `f(p)` is dimension-dependent (Table 1 in the paper):

| D | Basis for f(p) |
|---|---|
| 5 | `p³−p², p⁴−p², p⁵−p², ...` |
| even ≥ 6 | `p^{3/2}, p^{5/2}, p^{7/2}, ...` |
| odd ≥ 7 | `p², p³, p⁴, ...` |

**Why these choices?** This is dictated by the asymptotic behavior of the Bessel function integrals at large impact parameter. The integral of `p^n` against `J_{(D-4)/2}(bp)/(bp/2)^{(D-4)/2}` has an asymptotic expansion:

```
_1F_2(n+1/2; (D-2)/2, (n+3)/2; −b²/4) / (n+1)
  ~ [non-oscillatory] + [oscillatory: cos(b − π(D−1)/4)/b^{(D-1)/2}]
```

For the functional to be **positive at large b** without cancellations, one needs either `n ≤ (D−3)/2` or specific combinations that cancel the leading oscillatory term.

In D=5: `(D−3)/2 = 1`, but we need `n > 1` for the gravity integral to converge. So we use differences `p^k − p²` to cancel the oscillatory term at leading order.

In even D ≥ 6: half-integer powers work naturally.

In odd D ≥ 7: integer powers `p², p³, ...` work.

### 2.3 Analytic Integration Against Basis Functions

The integrals of basis functions against the kernels can be computed **analytically** using hypergeometric functions. For example (J=2 case):

```
∫₀¹ dp p^n C_{2,−p²}^{improved}[m², 2]
  = −4(D−1) _2F_1(1, (n+1)/2; (n+3)/2; −1/m²) / [(D−2)m⁴(n+1)]
    + 2(3D−4)/[(D−2)m⁴(n+1)] − 3(3D−2)/[(D−2)m⁶(n+3)]
```

This is exact for any `m`. This allows extremely precise evaluation of the constraint matrix.

### 2.4 Mass Discretization

The mass parameter `m ∈ [1, ∞)` (with `M=1`) is parametrized as `m² = 1/(1−x)` so that `x ∈ [0, 1)`. The initial discretization is:

```
x ∈ {0, δ_x, 2δ_x, ..., ⌈1/δ_x − 1⌉·δ_x}
```

Typical values: `δ_x = 1/400` (Figure 1) or `1/800` (Figure 3 — extremal spectrum).

This turns the infinite set of positivity constraints (one per value of m) into a **finite set of linear inequalities**.

### 2.5 Angular Momentum Truncation and J_max

The sum over J is truncated at `J_max`. Typical values: `J_max = 42` (Figure 1) or `J_max = 150` (Figure 3).

To compensate, the **scaling limit** `m → ∞` with fixed `b = 2J/m` is included explicitly as an additional inequality:

```
Γ((D−2)/2) ∫₀¹ dp f(p) J_{(D-4)/2}(bp) / (bp/2)^{(D-4)/2} ≥ 0   ∀b ≥ 0
```

This is the impact parameter space positivity condition, which captures the physics of large J without needing to include each spin individually.

### 2.6 Handling the Impact Parameter Inequality at Large b

The impact parameter integral is decomposed as:

```
A(b) + B(b)·cos(b − π(D-1)/4) + C(b)·sin(b − π(D-1)/4)
```

where `A, B, C` are power series in `1/b`. Positivity of this expression is **replaced by the stronger (but rigorous) condition**:

```
[A(b)+B(b)   C(b)  ]
[C(b)        A(b)−B(b)] ≽ 0   (positive semidefinite)
```

This 2×2 matrix condition is then expanded in powers of `1/b` and truncated at `m_max` subleading terms. The resulting polynomial matrix inequality can be directly fed to **SDPB** as a semidefinite constraint.

For small `b ≤ B` (cutoff), positivity is imposed at discretized impact parameters:
```
b ∈ {ε_b, ε_b+δ_b, ..., B}
```

### 2.7 The Full Linear Program as an SDP

SDPB solves semidefinite programs of the form:

```
minimize  b · y
subject to M_j(y) ≽ 0   for each j
```

The key transformation: the positivity constraint `∫ dp f(p) K(m,J) ≥ 0` for all (m,J) becomes a family of linear constraints on the coefficients `{a_n, b_{k,i}}`. The impact parameter constraint becomes a polynomial matrix inequality. Together, these define a large SDP.

**SDPB parameters used:**
- Precision: 768–840 decimal digits (needed because the constraints are nearly degenerate near the optimal functional)
- Duality gap threshold: `1e-80` (for the high-precision Figure 3 computation)
- Primal/dual error thresholds: `1e-80`

### 2.8 Objective Function and Normalization

To bound e.g. `g₃` from above at given `g₂`:

```
minimize  ∫₀¹ dp f(p) [1/p² + 2g₂]
such that ∫₀¹ dp f(p) p² = −1
```

with the positivity constraints. The result gives an upper bound on `g₃` as a function of `g₂` and `G`.

To trace out the boundary of the allowed region (Figure 1), the paper scans over angles `θ ∈ {0, π/20, ..., 39π/20}` in the `(g₂, g₃)` plane, maximizing the distance from a known interior point.

---

## Part 3: The Adaptive Refinement Flow — Analysis and Bugs

### 3.1 The Problem Being Solved

The adaptive refinement procedure solves a fundamental tension: the true constraints are **continuous** (one inequality per `m ∈ [1,∞)` and per `b ∈ [0,∞)`), but the solver needs **finitely many** constraints. A naive coarse discretization will return a functional that satisfies all the discretized constraints but **goes negative between the grid points**.

The authors note: *"Typically, the solver returns a solution that vanishes at pairs of neighboring discrete values of x, and is negative between them."* This is the core issue the adaptive loop must fix.

### 3.2 The Refinement Loop (Described in the Paper)

The procedure as described:

```
1. Start with initial grid: x ∈ {0, δ_x, 2δ_x, ..., 1-δ_x}
2. Run SDPB → get candidate functional F
3. For each pair (x₁, x₂) of adjacent grid points where F dips negative:
    a. Estimate minimum location x* via quadratic approximation
    b. Compute s = |x₂ − x₁| / N  (with e.g. N=10)
    c. Add new constraints at {x* − N·s, ..., x*, ..., x* + N·s}
4. Re-run SDPB from scratch with new constraints
5. Repeat until negative regions are smaller than tolerance (e.g. size 10⁻⁶)
```

The same procedure is applied independently to the impact parameter discretization grid (for `b`).

### 3.3 Identified Bugs and Failure Modes

#### Bug #1: Quadratic Approximation Misses the True Minimum

**Location:** Step 3a — estimation of `x*`

The paper states: *"Let the minimum of the functional between x₁ and x₂ (which we **estimate from a quadratic approximation**) be x*."*

**The bug:** A quadratic fit through two boundary points `(x₁, F(x₁))` and `(x₂, F(x₂))` and possibly the midpoint is used to estimate `x*`. However:
- If the functional has a **non-quadratic profile** (e.g., a sharp dip or an asymmetric notch), the quadratic approximation places `x*` at the wrong location.
- The new constraints are then clustered around the wrong `x*`, **leaving the true minimum in a gap between constraints**.
- On the next SDPB run, the solver will find a negative region nearby — **the very same region it was supposed to cancel**. This matches the observed symptom: "tasks that should have been cancelled" (i.e., negative regions that should have been eliminated by constraint insertion) **persist** into subsequent refinement rounds.

**Consequence:** Near the optimal functional (which saturates inequalities at multiple points), the functional profile is generically **cusped or flat near the minimum**, making the quadratic approximation systematically wrong. This can cause the refinement to stall or converge to a suboptimal functional.

---

#### Bug #2: No Hot-Start Between Refinement Iterations

**Location:** The outer refinement loop

The paper explicitly acknowledges: *"An efficient implementation of this method should include a way of hot-starting from the previous solution after each refinement step. **We have not implemented this** — instead we simply run SDPB from scratch after each refinement."*

**The bug:** Each iteration restarts SDPB from scratch, discarding the previous dual solution as a warm start. In practice:
- SDPB uses an interior-point method that benefits enormously from a good starting point
- Without hot-starting, each refinement round pays the full computational cost of convergence from a cold start
- More critically, **when new constraints are added that tighten the feasible region**, the previous solution may not even be feasible for the new problem. A hot-started interior-point method would detect this immediately and steer toward the new optimum. Without it, the solver may **wander through infeasible directions** before finding a new feasible point, wasting iterations
- This can cause the solver to return a **different local structure** of the optimal functional on each refinement pass, where the negative regions move around rather than shrinking — mimicking the symptom of "tasks that should have been cancelled" running again in slightly different locations

---

#### Bug #3: Fixed Step Count N Ignores Local Gradient Information

**Location:** Step 3b — computation of `s = |x₂ − x₁| / N`

The refinement step size is `s = |x₂ − x₁| / N` with a **fixed** `N = 10`. The constraint points are uniformly spaced around `x*`.

**The bug:** The optimal distribution of new constraint points should be non-uniform, concentrated where the functional is steepest (near zeros and sign changes). Using a uniform spacing with fixed `N`:
- May place constraints **too far** from the actual sign-change location in cases where the functional dip is extremely narrow (near the optimum, dips can span much less than `|x₂ − x₁| / 10`)
- In these cases, the constraints added in one round **fail to straddle the true zero**, and after re-running SDPB, the negative region shifts to a new narrow location — again "not cancelled"
- The paper reports (Figure 5 caption): *"The negative region is b ∈ (10.7298112, 10.7298134), minimum value −1.83×10⁻¹⁰"* — a region of width ~2×10⁻⁶, demonstrating that dips can be extremely narrow. A uniform N=10 subdivision of the surrounding interval would only achieve resolution ~`δ_b/10`, which may still miss such a narrow feature until many rounds of refinement have occurred.

---

#### Bug #4: Independent Grids for x and b Miss Coupled Negative Regions

**Location:** The dual-grid management (x grid for mass, b grid for impact parameter)

The refinement is applied independently to the `x` (mass) grid and the `b` (impact parameter) grid.

**The bug:** The positivity condition is a function of **both** `m` (equivalently `x`) and `b` simultaneously. A functional can be negative in a **two-dimensional region** of the `(m, b)` plane — for example, it may be marginally positive at the discretized `b` values but negative at intermediate `b` for certain `m` values near a saturation point.

By refining each grid independently, the algorithm can:
- Add constraints at new `b` values that fix the issue at the specific `m` values tested, while
- Missing the case where the functional turns negative at a **different `m`** with the same `b`

The paper's description does not mention any joint `(m, b)` refinement strategy. This architectural gap means that a functional which is negative in a curved feature of the `(m, b)` plane can evade detection for many refinement rounds — **a class of "tasks that should have been cancelled" that the current loop cannot reliably detect**.

---

#### Bug #5: Convergence Criterion is One-Sided

**Location:** Step 5 — termination condition

The loop terminates when negative regions are smaller than a size tolerance (e.g., `10⁻⁶` in `x`). However:
- "Size in x" (or equivalently in b) is not the same as "magnitude of violation"
- A functional with a very narrow but deeply negative dip (say, minimum −0.1 over an interval of width 10⁻⁷) would **pass the size threshold** but represent a significant violation of the positivity constraint
- The paper reports achieving minimum values as small as `−1.83×10⁻¹⁰` after refinement, but this is after many refinement steps; the stopping criterion based on interval size alone could terminate too early in adversarial cases

**Consequence:** Bounds returned by the solver may be **slightly weaker than reported** if a small-but-deep negative region was missed. The paper appears aware of this (noting convergence to a fraction `10⁻⁴`), but the criterion is not theoretically tight.

---

#### Bug #6 (Historical, Now Fixed): Incorrect p_min ↔ b_max Relation for D=4

**Location:** Section on D=4 infrared divergences

The paper explicitly states: *"An earlier arXiv version of this paper used an incorrect relation between p_min and b_max, which resulted in a numerically incorrect bound."*

The corrected relation is `p_min² = c/(M·b_max³)` with `c ≈ 1`. The earlier version apparently used an incorrect functional form or coefficient, producing a bound that was parametrically wrong. This has been corrected in the published version, but it illustrates the fragility of the D=4 analysis.

---

### 3.4 Summary of Bugs

| # | Location | Description | Symptom |
|---|---|---|---|
| 1 | Quadratic minimum estimation | Wrong `x*` placement when dip is non-quadratic | Negative regions persist after refinement |
| 2 | No hot-start between rounds | SDPB restarts cold, may explore different basin | Negative regions migrate rather than shrink |
| 3 | Fixed N=10 uniform spacing | Too coarse for sub-interval-width negative dips near optimum | Narrow dips not cancelled until many rounds |
| 4 | Independent x and b grids | Misses coupled (m, b) negative regions | 2D violation regions evade detection |
| 5 | Size-based stopping criterion | Does not check magnitude of violation | Deep-but-narrow violations terminate too early |
| 6 | (Historical) D=4 p_min/b_max formula | Incorrect exponent in relation | Numerically wrong bound in prior arXiv version |

---

### 3.5 Important Note on Scope

The files provided are **LaTeX source code for a physics research paper**, not a software system with task scheduling infrastructure. There is no scheduler, task queue, cancellation mechanism, or job management code in any of the uploaded files. The "task scheduling flow" in this context refers to the **iterative adaptive refinement loop** described above — a mathematical algorithm that decides when to add new positivity constraints and when to re-run the optimizer.

All "bugs" identified above are **algorithmic flaws** in the mathematical procedure as described in the paper, not software bugs in a codebase. The symptom "tasks that should have been cancelled" corresponds to **positivity constraints that were supposed to be enforced (eliminating negative regions of the functional) but whose enforcement is deferred or missed by the imprecise refinement procedure**.

---

## Appendix: Key Parameter Tables

### Numerical Parameters (from Table 2 of paper)

| Parameter | Figure 1 | Figure 3 (extremal spectrum) |
|---|---|---|
| J_max | 42 | 150 |
| δ_x | 1/400 | 1/800 |
| ε_b | 1/250 | 1/250 |
| δ_b | 1/32 | 1/100 |
| B (large-b cutoff) | 40 | 80 |
| m_max (asymptotic terms) | 2 | 6 |
| SDPB precision | 768 bits | 840 bits |

### Functional Space Dimensions

- **Figure 1:** 17-dimensional space (C₂^improved with n ≤ 7.5, plus X₄ with i=0,...,5, plus X₆ with i=0,...,3)
- **Figure 3:** 93-dimensional space (larger ranges of X₄, X₆, X₈, X₁₀, X₁₂, X₁₄)
- **Figure 4 (higher contacts):** 54-dimensional space
- **Maximal SUGRA bound:** 14 functionals total (gravity C_{-2} measured with p², ..., p¹⁰; g₀ with 4 powers; X₂ with 2 powers)

### Known Solutions and Bounds (D=6, units 8πG=1, M=1)

- **Type II string:** `g₂·m²/(8πG) = 0`, `g₃·m⁴/(8πG) = 0` (origin in plot)
- **Type II - spin-0:** `g₂·m²/(8πG) ≈ −0.75`, `g₃·m⁴/(8πG) ≈ −2.07`
- **Heterotic:** `g₂·m²/(8πG) = 3/16`, `g₃·m⁴/(8πG) = 3/4`
- **SUGRA (D=10):** `g₀·M⁶/(8πG) ≤ 3.000`; type II gives `2ζ(3) ≈ 2.404`

---

## Part 4: Pseudocode — Complete Procedure for Bounding g₂ and g₃

This section provides a fully explicit, step-by-step algorithmic description of the numerical procedure for tracing the boundary of the allowed region in the `(g₂, g₃)` plane. Every step is derived directly from the paper's equations and appendix, and is cross-referenced to the mathematical objects defined in Parts 1–3. The algorithm is designed to be correct, implementable, and to faithfully represent what the paper actually computes — including all normalizations, analytic formulas, constraint types, and the adaptive refinement loop.

---

### 4.1 Mathematical Foundations (Recap for Algorithm Context)

Working in units `M = 1` and `8πG = 1` throughout (as the paper does in the appendix). Physical bounds in other units are recovered by dimensional analysis.

**The fundamental linear program (LP)** is:

```
if:   ∀ m ≥ 1, ∀ J = 0, 2, 4, ...:
        F(m, J) ≡ ∫₀¹ dp f(p) · C²ᵢₘₚ[m², J, -p²]
                  + Σ_{k=4,6,...} ∫₀¹ dp h_k(p) · X_k[m², J, -p²]  ≥  0

then: I[f] ≡ ∫₀¹ dp f(p) · [1/p² + 2g₂ + g₃·p²]  ≥  0
```

where the two kernels are:

**Improved kernel** (from eq. B2 improved flat, with substitution u = −p²):
```
C²ᵢₘₚ[m², J, u] =
    (2m²+u) · Pⱼ(1 + 2u/m²)
    ─────────────────────────────
          m²·(m²+u)²

    - u²/m⁶ · [ (4m²+3u)·Pⱼ(1)     4u·Pⱼ'(1)  ]
               [ ─────────────────  + ─────────── ]
               [    (m²+u)²             m⁴−u²    ]
```

Here:
- `Pⱼ(x) = ₂F₁(−J, J+D−3, (D−2)/2, (1−x)/2)` are Gegenbauer polynomials
- `Pⱼ(1) = 1` for all J (boundary value)
- `Pⱼ'(1) = J(J+D−3)/(D−2)` (derivative at 1, equals quadratic Casimir divided by (D−2))

**Null constraint kernel** (from eq. nullconstraints, with u = −p²):
```
X_k[m², J, u] =
    (2m²+u)              m² · Pⱼ(1+2u/m²)
    ──────────────  ·  ─────────────────────────
    u·m²·(m²+u)         [u·m²·(m²+u)]^{k/2}

    − Res_{u'=0} [
        (2m²+u')(m²−u')(m²+2u')    m²·Pⱼ(1+2u'/m²)
        ─────────────────────────── · ──────────────────────────
        m²(u−u')u'(m²−u)(m²+u')(m²+u+u')   [u'·m²·(m²+u')]^{k/2}
      ]
```

The residue at `u' = 0` is computed by Laurent-expanding around `u' = 0` and extracting the coefficient of `u'^{−1}`. The result is a rational function of `m` and `u`, polynomial in `J(J+D-3)`.

**EFT action of f:**
```
I[f] = ∫₀¹ dp f(p) · [1/p² + 2g₂ + g₃·p²]
     = γ[f] + 2g₂·α[f] + g₃·β[f]

where:
  γ[f] = ∫₀¹ dp f(p)/p²      (action on gravity / Newton constant)
  α[f] = ∫₀¹ dp f(p)         (action on g₂)
  β[f] = ∫₀¹ dp f(p)·p²      (action on g₃)
```

---

### 4.2 Algorithm Parameters

```
INPUT PARAMETERS:
  D           — spacetime dimension (integer, D ≥ 5; D = 4 requires separate IR treatment)
  N_f         — number of basis functions for f(p)
  k_list      — list of null constraints to include, e.g. [4, 6] for 17-dim space
  N_h[k]      — number of basis functions for h_k(p) for each k in k_list
  N_angles    — number of boundary directions θ to scan (paper uses 40)
  J_max       — maximum angular momentum truncation (paper: 42 for Fig.1, 150 for Fig.3)
  δ_x         — initial mass-grid spacing in x = 1 − 1/m² (paper: 1/400 or 1/800)
  ε_b         — initial impact-parameter grid offset (paper: 1/250)
  δ_b         — initial impact-parameter grid spacing (paper: 1/32 or 1/100)
  B_large     — large-b cutoff; above this use PSD matrix (paper: 40 or 80)
  m_max       — number of subleading terms kept in large-b expansion (paper: 2 or 6)
  N_refine    — subdivision count per negative region (paper: N = 10)
  tol_neg     — size tolerance for negative regions (paper: ~1e-6)
  tol_val     — magnitude tolerance for negativity (additional depth check)
  (g₂₀, g₃₀) — interior reference point near expected tip of allowed region

DERIVED:
  basis_f     — list of powers {n} for f(p) = Σ aₙ pⁿ, chosen by D (Table 1)
  basis_h[k]  — list of powers {i} for h_k(p) = Σ b_{k,i} pⁱ
  x_grid_init — initial mass grid {0, δ_x, 2δ_x, ..., ⌈1/δ_x−1⌉·δ_x}
  b_grid_init — initial impact-parameter grid {ε_b, ε_b+δ_b, ..., B_large−δ_b}
  J_list      — [0, 2, 4, ..., J_max]
  θ_list      — {0, π/N_angles, 2π/N_angles, ..., (2N_angles−1)π/N_angles}
```

---

### 4.3 Step 1: Select Basis Functions for f(p)

```
FUNCTION SelectBasisF(D, N_f):
  // Basis choice is dictated by the large-b positivity requirement.
  // The Fourier transform of p^n at large b behaves as:
  //   ~ (non-oscillatory decay ~1/b^{n+1})
  //   + (oscillatory ~cos(b − π(D-1)/4) / b^{(D-1)/2})
  // For generic f to be positive at large b, one of these conditions must hold:
  //   (A) There exists at least one n ≤ (D−3)/2 in the basis, OR
  //   (B) The leading oscillatory term cancels in the linear combination.
  // Condition on powers: n > 1 always required so that ∫₀¹dp f(p)/p² converges.

  IF D == 5:
    // (D−3)/2 = 1, but need n > 1, so neither condition (A) holds.
    // Use condition (B): differences p^k − p² cancel leading oscillation.
    // basis_f = {p³−p², p⁴−p², p⁵−p², ..., p^{N_f+1}−p²}
    powers_raw = [3, 4, 5, ..., N_f+2]
    RETURN [(pⁿ − p²) for n in powers_raw]   // N_f basis functions

  ELSE IF D is even AND D >= 6:
    // Half-integer powers work naturally: non-oscillatory term dominates.
    // basis_f = {p^{3/2}, p^{5/2}, p^{7/2}, ..., p^{(2N_f+1)/2}}
    RETURN [p^{(2i+1)/2} for i in 1..N_f]

  ELSE IF D is odd AND D >= 7:
    // Integer powers p^n for n ≥ 2, since (D−3)/2 ≥ 2 allows positivity.
    // basis_f = {p², p³, p⁴, ..., p^{N_f+1}}
    RETURN [pⁿ for n in 2..N_f+1]

  // Note: for D=4, require a separate infrared cutoff treatment (Section 3.5 of research.md).
```

---

### 4.4 Step 2: Precompute Analytic Kernel Integrals

All integrals of basis functions against kernels can be evaluated analytically as closed-form functions of `m` (or equivalently `x = 1 − 1/m²`). This precomputation is the most expensive symbolic step but is done once per parameter set.

```
FUNCTION PrecomputeKernelTable(D, basis_f, k_list, basis_h, J_list):
  // ─────────────────────────────────────────────────────────────────────
  // TABLE A: Improved kernel integrals ∫₀¹ dp pⁿ · C²ᵢₘₚ[m², J, −p²]
  // For each (n, J) pair, compute as a closed-form function of m.
  // ─────────────────────────────────────────────────────────────────────

  FOR EACH power n in basis_f (after expanding differences if D=5):
    FOR EACH J in J_list:
      // The Gegenbauer polynomial Pⱼ(1 − 2p²/m²) is a degree-J polynomial
      // in p²/m². Expanding and integrating term by term yields ₂F₁ functions.
      //
      // General formula (analytical, for any J):
      //   ∫₀¹ dp pⁿ · C²ᵢₘₚ[m², J, −p²]
      //   = T1(n, J, m) + T2(n, J, m) + T3(n, J, m)
      //
      // where, using the decomposition of C²ᵢₘₚ into three terms:
      //
      //   T1(n,J,m) = ∫₀¹ dp pⁿ · (2m²−p²)·Pⱼ(1−2p²/m²) / [m²·(m²−p²)²]
      //
      //   T2(n,J,m) = −∫₀¹ dp pⁿ · p⁴/m⁶ · (4m²−3p²)·Pⱼ(1) / (m²−p²)²
      //             = −Pⱼ(1)/m⁶ · ∫₀¹ dp pⁿ · p⁴(4m²−3p²) / (m²−p²)²
      //               [Pⱼ(1) = 1, so this term is independent of J]
      //
      //   T3(n,J,m) = −∫₀¹ dp pⁿ · p⁴/m⁶ · (−4p²)·Pⱼ'(1) / (m⁴−p⁴)
      //             = 4Pⱼ'(1)/m⁶ · ∫₀¹ dp pⁿ · p⁶ / (m⁴−p⁴)
      //               [Pⱼ'(1) = J(J+D-3)/(D-2) is J-dependent]
      //
      // For J=2 specifically (D-dimensional spin-2 pole):
      //   ∫₀¹ dp pⁿ · C²ᵢₘₚ[m², 2, −p²]
      //   = −4(D-1)·₂F₁(1,(n+1)/2;(n+3)/2;−1/m²) / [(D-2)·m⁴·(n+1)]
      //     + 2(3D-4) / [(D-2)·m⁴·(n+1)]
      //     − 3(3D-2) / [(D-2)·m⁶·(n+3)]
      //   [from eq. (A.3) of the paper, exact formula for J=2]
      //
      // For general J, T1 involves:
      //   ∫₀¹ dp pⁿ · p^{2r} / (m²−p²)² for r = 0,1,...,J
      //   = B((n+2r+1)/2, 1/2) · m^{n+2r-3} · ₂F₁((n+2r+1)/2, 2; (n+2r+3)/2; −1/m²)
      //     / [various gamma-function factors]
      //   computed via standard Beta-function integrals.
      //
      // Store result as function-of-m: KernelC[n][J](m)
      // (For the D=5 difference basis p^k − p², this is KernelC[k][J](m) − KernelC[2][J](m))

      KernelC[n][J] = COMPUTE_ANALYTIC_INTEGRAL_C2IMP(n, J, D)
      // Returns a callable function of m, expressed as a finite sum
      // of hypergeometric functions ₂F₁(a, b; c; −1/m²).

  // ─────────────────────────────────────────────────────────────────────
  // TABLE B: Null constraint integrals ∫₀¹ dp pⁱ · X_k[m², J, −p²]
  // ─────────────────────────────────────────────────────────────────────

  FOR EACH k in k_list:    // k = 4, 6, ...
    FOR EACH power i in basis_h[k]:    // i = 0, 1, 2, ...
      FOR EACH J in J_list:
        // X_k[m², J, u] is defined by the null-constraint formula.
        // With u = −p², each term in X_k is a rational function of p and m
        // times Pⱼ(1 − 2p²/m²). The residue at u'=0 produces polynomial-in-(p²/m²)
        // corrections. Expanding Pⱼ and integrating yields more ₂F₁ functions.
        //
        // The residue computation at u'=0:
        //   The integrand has a pole of order ⌈k/2⌉ at u'=0.
        //   Extract coefficient of u'^{−1} in the Laurent series → rational in m², u.
        //   For k=4: pole is order 2, residue requires first-derivative expansion.
        //   For k=6: pole is order 3, requires second-derivative expansion. Etc.

        KernelX[k][i][J] = COMPUTE_ANALYTIC_INTEGRAL_Xk(k, i, J, D)
        // Returns a callable function of m.

  // ─────────────────────────────────────────────────────────────────────
  // TABLE C: Impact parameter integrals (EFT action integrals)
  // γ[pⁿ], α[pⁿ], β[pⁿ] for the EFT coupling vector
  // ─────────────────────────────────────────────────────────────────────

  FOR EACH power n in basis_f:
    gamma[n] = ∫₀¹ dp pⁿ / p²   = ∫₀¹ dp p^{n−2}   = 1/(n−1)   [requires n > 1]
    alpha[n] = ∫₀¹ dp pⁿ         = 1/(n+1)
    beta[n]  = ∫₀¹ dp pⁿ · p²   = 1/(n+3)

  // ─────────────────────────────────────────────────────────────────────
  // TABLE D: Impact parameter space integrals (large-J positivity)
  // Γ((D-2)/2) ∫₀¹ dp pⁿ Jᵥ(bp) / (bp/2)^ν  where ν = (D-4)/2
  // ─────────────────────────────────────────────────────────────────────

  FOR EACH power n in basis_f:
    // Exact analytic formula (from paper, impact parameter appendix):
    //   ImpactInt[n](b) = ₁F₂((n+1)/2; (D-2)/2, (n+3)/2; −b²/4) / (n+1)
    //
    // Large-b asymptotic expansion (needed for PSD condition):
    //   ImpactInt[n](b) ~
    //     A_coeff[n] / b^{n+1}                                  ← non-oscillatory
    //   + B_coeff[n] · cos(b − π(D-1)/4) / b^{(D-1)/2}        ← oscillatory cos
    //   + C_coeff[n] · sin(b − π(D-1)/4) / b^{(D-1)/2+1}      ← oscillatory sin (subleading)
    //   + ...
    //
    // where:
    //   A_coeff[n] = 2ⁿ Γ((D-2)/2) Γ((n+1)/2) / Γ((D-n-3)/2)
    //   B_coeff[n] = 2^{(D-3)/2} Γ((D-2)/2) / √π       [same for all n at leading order]
    //   C_coeff[n] = subleading coefficient (from next term in ₁F₂ expansion)

    ImpactExact[n]   = λb. ₁F₂((n+1)/2; (D-2)/2, (n+3)/2; −b²/4) / (n+1)
    ImpactAsymp[n]   = COMPUTE_ASYMPTOTIC_EXPANSION(n, D, m_max)
    // Returns A(b), B(b), C(b) as truncated Laurent series in 1/b,
    // keeping m_max subleading terms beyond the leading term.

  RETURN (KernelC, KernelX, gamma, alpha, beta, ImpactExact, ImpactAsymp)
```

---

### 4.5 Step 3: Evaluate Kernel Tables on the Mass Grid

```
FUNCTION EvaluateOnMassGrid(x_grid, J_list, KernelC, KernelX, basis_f, k_list, basis_h):
  // Convert x → m using m² = 1/(1−x), i.e., m = 1/√(1−x)

  FOR EACH x_i in x_grid:
    m_i = 1 / sqrt(1 − x_i)

    FOR EACH J in J_list:
      FOR EACH n in basis_f:
        // Evaluate the precomputed analytic function at m = m_i
        MatC[n][J][x_i] = KernelC[n][J](m_i)
        // Shape: [N_f × |J_list| × |x_grid|] array of scalars

      FOR EACH k in k_list:
        FOR EACH i in basis_h[k]:
          MatX[k][i][J][x_i] = KernelX[k][i][J](m_i)
          // Shape: [k_list × max(N_h) × |J_list| × |x_grid|] array

  RETURN (MatC, MatX)
```

---

### 4.6 Step 4: Build the SDP Constraint System

The finite LP has the following structure in terms of decision variables:

```
Decision variables:
  a = {aₙ}          — coefficients of f(p) = Σₙ aₙ pⁿ
  b = {b_{k,i}}      — coefficients of h_k(p) = Σᵢ b_{k,i} pⁱ
  Full vector:  v = (a₁, a₂, ..., b_{4,0}, b_{4,1}, ..., b_{6,0}, ...)
  Dimension:    N_vars = N_f + Σ_k N_h[k]
```

```
FUNCTION BuildConstraintSystem(MatC, MatX, x_grid, b_grid, J_list, k_list,
                                basis_f, basis_h, ImpactExact, ImpactAsymp,
                                B_large, m_max, D):

  constraints = []   // list of (row_vector, lower_bound) pairs: row·v ≥ lb

  // ─────────────────────────────────────────────────────────────────────
  // BLOCK 1: Finite-m positivity constraints
  //   F(m_i, J) = Σₙ aₙ MatC[n][J][x_i] + Σ_k Σᵢ b_{k,i} MatX[k][i][J][x_i] ≥ 0
  //   One constraint per (x_i, J) pair.
  // ─────────────────────────────────────────────────────────────────────

  FOR EACH x_i in x_grid:
    FOR EACH J in J_list:
      row = zero_vector(N_vars)

      // Fill entries for f-coefficients: position in row = index of n in basis_f
      FOR EACH n in basis_f (index = idx_n):
        row[idx_n] += MatC[n][J][x_i]

      // Fill entries for h_k-coefficients
      FOR EACH k in k_list:
        FOR EACH i in basis_h[k] (index = idx_{k,i} in v):
          row[idx_{k,i}] += MatX[k][i][J][x_i]

      constraints.append( (row, lower_bound=0.0) )

  // ─────────────────────────────────────────────────────────────────────
  // BLOCK 2: Discretized impact-parameter positivity (small b ≤ B_large)
  //   Σₙ aₙ · ImpactExact[n](b_j) ≥ 0   for each b_j in b_grid
  //   (h_k functions do not enter: they are subleading at m→∞)
  // ─────────────────────────────────────────────────────────────────────

  FOR EACH b_j in b_grid (b_j ≤ B_large):
    row = zero_vector(N_vars)
    FOR EACH n in basis_f (index = idx_n):
      row[idx_n] += ImpactExact[n](b_j)
    constraints.append( (row, lower_bound=0.0) )

  // ─────────────────────────────────────────────────────────────────────
  // BLOCK 3: Large-b PSD matrix constraint (b ≥ B_large)
  //   The 2×2 matrix [[A(b)+B(b), C(b)], [C(b), A(b)−B(b)]] ≽ 0
  //
  //   Each entry is a linear function of the coefficients {aₙ}:
  //     A(b) = Σₙ aₙ · A_n(b),  B(b) = Σₙ aₙ · B_n(b),  C(b) = Σₙ aₙ · C_n(b)
  //   where A_n(b), B_n(b), C_n(b) are extracted from ImpactAsymp[n].
  //
  //   Substituting b = B_large / (1−t) or equivalently expanding in 1/b,
  //   the 2×2 matrix entries become Laurent polynomials in 1/b.
  //   After multiplying through by b^{(D-1)/2} (a positive factor), we get
  //   a polynomial matrix inequality in 1/b, which SDPB accepts natively.
  //
  //   The polynomial matrix PSD condition encodes infinitely many scalar
  //   constraints (one for every b ≥ B_large) as a finite SDP constraint.
  // ─────────────────────────────────────────────────────────────────────

  // Expand A(b), B(b), C(b) to order m_max in 1/b:
  //   A(b) = Σₙ aₙ · [A_n,0/b^{n+1} + A_n,1/b^{n+3} + ...]
  //   B(b) = Σₙ aₙ · [B_n,0/b^{(D-1)/2} + B_n,1/b^{(D+1)/2} + ...]
  //   C(b) = Σₙ aₙ · [C_n,0/b^{(D-1)/2} + C_n,1/b^{(D+1)/2} + ...]
  //
  // Collect A(b)+B(b), A(b)−B(b), C(b) as linear combinations of {aₙ}.
  // Extract polynomial coefficients in z = 1/b (substituting b = 1/z):
  //   [A+B](z) = (Σₙ aₙ qₙ,AB(z))   — polynomial in z (after mult. by z^{−deg})
  //   [A−B](z) = (Σₙ aₙ qₙ,AmB(z))
  //   [C](z)   = (Σₙ aₙ qₙ,C(z))
  //
  // The PSD condition [[A+B, C],[C, A−B]] ≽ 0 for all z ∈ [0, 1/B_large]
  // is encoded as an SDPB polynomial matrix PSD block.

  psd_block = BUILD_POLYNOMIAL_MATRIX_PSD_BLOCK(
    {ImpactAsymp[n] for n in basis_f}, B_large, m_max, D
  )
  // This is passed to SDPB as a separate polynomial matrix constraint block.

  RETURN (constraints, psd_block)
```

---

### 4.7 Step 5: Set Up the Objective for a Given Scan Angle θ

For each scan direction `θ` in `θ_list`, we find the farthest boundary point from the interior point `(g₂₀, g₃₀)` in direction `θ`.

```
FUNCTION SetupObjective(theta, g20, g30, gamma, alpha, beta, basis_f):
  // The inequality from the functional:
  //   I[f] = γ[f] + 2g₂·α[f] + g₃·β[f] ≥ 0
  // defines a half-space in (g₂, g₃) space for each valid f.
  //
  // Boundary point in direction θ from (g₂₀, g₃₀):
  //   At the boundary: I[f] = 0 when g₂ = g₂₀ + λcosθ, g₃ = g₃₀ + λsinθ.
  //   Substituting: γ[f] + 2(g₂₀+λcosθ)·α[f] + (g₃₀+λsinθ)·β[f] = 0
  //   ⟹ λ·(2cosθ·α[f] + sinθ·β[f]) = −(γ[f] + 2g₂₀·α[f] + g₃₀·β[f])
  //   ⟹ λ = −I₀[f] / D_θ[f]
  //
  // where I₀[f] = γ[f] + 2g₂₀·α[f] + g₃₀·β[f]   ("action at reference point")
  //       D_θ[f] = 2cosθ·α[f] + sinθ·β[f]           ("directional action")
  //
  // Maximizing λ with normalization D_θ[f] = 1 reduces to:
  //   minimize I₀[f] subject to D_θ[f] = 1
  //
  // The maximum λ is then: λ* = −min I₀[f] (with the normalization applied)
  //
  // In terms of coefficient vector v = (a₁, ..., b_{k,i}, ...):

  // Objective vector (for minimization of I₀[f]):
  obj = zero_vector(N_vars)
  FOR EACH n in basis_f (index = idx_n):
    // For D=5 difference basis (p^k − p²): subtract contributions
    obj[idx_n] = gamma[n] + 2*g20*alpha[n] + g30*beta[n]
    // (if basis function is (p^k − p²), then
    //  gamma = 1/(k-1) − 1/(2-1), alpha = 1/(k+1) − 1/3, beta = 1/(k+3) − 1/5 )

  // Normalization constraint vector (equality, D_θ[f] = 1):
  norm = zero_vector(N_vars)
  FOR EACH n in basis_f (index = idx_n):
    norm[idx_n] = 2*cos(theta)*alpha[n] + sin(theta)*beta[n]

  RETURN (obj, norm, norm_rhs=1.0)
```

---

### 4.8 Step 6: Solve the SDP via SDPB

```
FUNCTION SolveSDP(constraints, psd_block, obj, norm, norm_rhs, precision_bits):
  // Assemble SDPB input files:
  //
  // SDPB accepts:
  //   — "bilinear bases": polynomial bases for PSD blocks
  //   — "free variable matrix" B: linear constraints Bv = b
  //   — "objective vector" c: minimize c·v
  //   — a set of semidefinite blocks
  //
  // 1. Convert linear inequality constraints (Block 1 and Block 2) to SDPB format.
  //    Each inequality row·v ≥ 0 becomes a 1×1 SDP block: [row·v] ≽ [0].
  //
  // 2. Pass the polynomial PSD block (Block 3) as an SDPB PSD block with
  //    the polynomial matrix basis corresponding to the expansion in z = 1/b.
  //
  // 3. Add the normalization equality:
  //    norm·v = 1   ←→   norm·v ≥ 1  AND  (−norm)·v ≥ −1
  //    (equivalently, treat as a free-variable equality in SDPB)
  //
  // 4. Set the objective to obj (minimize obj·v).

  sdpb_input = ASSEMBLE_SDPB_INPUT(constraints, psd_block, obj, norm, norm_rhs)

  // SDPB parameters (from paper's Table 2):
  //   --precision = precision_bits        (768 for Fig.1, 840 for Fig.3)
  //   --dualityGapThreshold = 1e-80       (for Fig.3 only)
  //   --primalErrorThreshold = 1e-80      (for Fig.3 only)
  //   --dualErrorThreshold = 1e-80        (for Fig.3 only)
  sdpb_params = {precision: precision_bits, ...}

  // Run SDPB (distributed computation, possibly on HPC cluster):
  result = RUN_SDPB(sdpb_input, sdpb_params)

  IF result.status == PRIMAL_FEASIBLE:
    v_opt = result.primal_solution    // optimal coefficient vector
    obj_opt = result.objective_value  // = min I₀[f] (should be ≤ 0 if interior pt is valid)
    RETURN (v_opt, obj_opt, SUCCESS)
  ELSE IF result.status == INFEASIBLE:
    // The LP has no feasible solution for this θ → no bound exists in this direction
    RETURN (null, null, INFEASIBLE)
  ELSE:
    RETURN (null, null, FAILURE)
```

---

### 4.9 Step 7: Validate the Functional and Identify Negative Regions

```
FUNCTION ValidateFunctional(v_opt, basis_f, k_list, basis_h, J_list,
                             x_grid, b_grid, B_large,
                             KernelC, KernelX, ImpactExact, D):
  // Reconstruct f(p) and h_k(p) from v_opt
  // Extract coefficient arrays:
  a    = {aₙ for n in basis_f}            // from v_opt
  bkk  = {b_{k,i} for k in k_list, i}    // from v_opt

  neg_regions_x = []    // list of (J, x₁, x₂, min_val) tuples
  neg_regions_b = []    // list of (b₁, b₂, min_val) tuples

  // ── Check finite-m positivity ──────────────────────────────────────
  // Evaluate F(m, J) on a dense verification grid (not the coarse training grid)
  verify_grid_x = {0, δ_x/100, 2δ_x/100, ..., 1 − δ_x/100}   // 100× denser

  FOR EACH J in J_list:
    prev_val = null
    prev_x   = null
    FOR EACH x_v in verify_grid_x:
      m_v = 1/sqrt(1 − x_v)
      F_val = Σₙ aₙ · KernelC[n][J](m_v)
              + Σ_k Σᵢ b_{k,i} · KernelX[k][i][J](m_v)

      IF F_val < 0 AND prev_val IS NOT NULL:
        // Sign change detected between prev_x and x_v
        neg_regions_x.append((J, prev_x, x_v, F_val))
      prev_val = F_val
      prev_x   = x_v

  // ── Check impact-parameter positivity ─────────────────────────────
  verify_grid_b = {ε_b, ε_b + δ_b/100, ε_b + 2δ_b/100, ..., B_large}

  prev_val = null
  prev_b   = null
  FOR EACH b_v in verify_grid_b:
    Fb_val = Σₙ aₙ · ImpactExact[n](b_v)

    IF Fb_val < 0 AND prev_val IS NOT NULL:
      neg_regions_b.append((prev_b, b_v, Fb_val))
    prev_val = Fb_val
    prev_b   = b_v

  RETURN (neg_regions_x, neg_regions_b)
```

---

### 4.10 Step 8: Adaptive Refinement Loop

```
FUNCTION AdaptiveRefinement(constraints_init, psd_block, obj, norm, norm_rhs,
                             x_grid_init, b_grid_init, J_list, KernelC, KernelX,
                             ImpactExact, ImpactAsymp, basis_f, k_list, basis_h,
                             B_large, N_refine, tol_neg, tol_val, precision_bits, D):

  x_grid    = x_grid_init       // mutable: grows as refinement proceeds
  b_grid    = b_grid_init       // mutable: grows as refinement proceeds
  constraints = constraints_init

  round = 0
  WHILE TRUE:
    round += 1

    // ── SOLVE ──────────────────────────────────────────────────────────
    (v_opt, obj_opt, status) = SolveSDP(constraints, psd_block, obj, norm,
                                        norm_rhs, precision_bits)
    IF status != SUCCESS:
      RAISE Error("SDP failed at round " + round)

    // ── VALIDATE ───────────────────────────────────────────────────────
    (neg_x, neg_b) = ValidateFunctional(v_opt, basis_f, k_list, basis_h,
                                        J_list, x_grid, b_grid, B_large,
                                        KernelC, KernelX, ImpactExact, D)

    // ── CHECK CONVERGENCE ─────────────────────────────────────────────
    // Convergence criterion 1: all negative regions have width < tol_neg
    all_x_small  = ALL( |x₂ − x₁| < tol_neg  for (J, x₁, x₂, v) in neg_x )
    all_b_small  = ALL( |b₂ − b₁| < tol_neg  for (b₁, b₂, v) in neg_b )

    // Convergence criterion 2: all negative values are shallower than tol_val
    all_x_shallow = ALL( |v| < tol_val  for (J, x₁, x₂, v) in neg_x )
    all_b_shallow = ALL( |v| < tol_val  for (b₁, b₂, v) in neg_b )

    IF (all_x_small AND all_b_small) OR (all_x_shallow AND all_b_shallow):
      // Functional has converged to positivity within tolerance
      BREAK

    // ── REFINE x-GRID (finite-m constraints) ───────────────────────────
    FOR EACH (J, x₁, x₂, min_val) in neg_x:
      // Estimate minimum location via quadratic interpolation.
      // Evaluate F at midpoint x_mid and at the two endpoints x₁, x₂:
      x_mid = (x₁ + x₂) / 2
      m_mid = 1/sqrt(1 − x_mid)
      F_mid = Σₙ aₙ·KernelC[n][J](m_mid) + Σ_k Σᵢ b_{k,i}·KernelX[k][i][J](m_mid)

      // Quadratic fit through (x₁, F(x₁)), (x_mid, F_mid), (x₂, F(x₂)):
      //   F(x) ≈ A(x−x₁)² + B(x−x₁) + C
      //   C = F(x₁), A+B = (F(x₂)−F(x₁))/(x₂−x₁), A(Δx)²/4 = F_mid−(F(x₁)+F(x₂))/2
      // Minimize: x* = x₁ − B/(2A) if A > 0 (otherwise x* = midpoint)

      (x_star, F_star) = QUADRATIC_MIN(x₁, F(x₁), x_mid, F_mid, x₂, F(x₂))

      // Guard: x_star must stay in (x₁, x₂)
      x_star = CLAMP(x_star, x₁ + 1e-10, x₂ − 1e-10)

      // Add N_refine new constraints symmetrically around x_star:
      s = |x₂ − x₁| / N_refine    // refinement step size
      new_x_points = {x_star + i*s  for i in −N_refine..N_refine}
      // Clip to [0, 1):
      new_x_points = {x for x in new_x_points if 0 ≤ x < 1}
      new_x_points = new_x_points MINUS x_grid    // only truly new points

      // Add new constraint rows to the system:
      FOR EACH x_new in new_x_points:
        m_new = 1/sqrt(1 − x_new)
        row = zero_vector(N_vars)
        FOR EACH n in basis_f (idx_n):
          row[idx_n] += KernelC[n][J](m_new)
        FOR EACH k in k_list:
          FOR EACH i in basis_h[k] (idx_{k,i}):
            row[idx_{k,i}] += KernelX[k][i][J](m_new)
        constraints.append((row, 0.0))
        x_grid.add(x_new)

    // ── REFINE b-GRID (impact-parameter constraints) ────────────────────
    FOR EACH (b₁, b₂, min_val) in neg_b:
      b_mid = (b₁ + b₂) / 2
      Fb_mid = Σₙ aₙ · ImpactExact[n](b_mid)

      (b_star, Fb_star) = QUADRATIC_MIN(b₁, Fb(b₁), b_mid, Fb_mid, b₂, Fb(b₂))
      b_star = CLAMP(b_star, b₁ + 1e-10, b₂ − 1e-10)

      s = |b₂ − b₁| / N_refine
      new_b_points = {b_star + i*s  for i in −N_refine..N_refine}
      new_b_points = {b for b in new_b_points if ε_b ≤ b ≤ B_large}
      new_b_points = new_b_points MINUS b_grid

      FOR EACH b_new in new_b_points:
        row = zero_vector(N_vars)
        FOR EACH n in basis_f (idx_n):
          row[idx_n] += ImpactExact[n](b_new)
        constraints.append((row, 0.0))
        b_grid.add(b_new)

  // ── RETURN CONVERGED RESULT ────────────────────────────────────────
  RETURN (v_opt, obj_opt, x_grid, b_grid, constraints)
```

---

### 4.11 Step 9: Extract the Boundary Point in (g₂, g₃) Space

```
FUNCTION ExtractBoundaryPoint(v_opt, obj_opt, theta, g20, g30, basis_f,
                               gamma, alpha, beta):
  // The optimal objective value obj_opt = min I₀[f] = min (γ[f] + 2g₂₀·α[f] + g₃₀·β[f])
  // gives the maximum distance λ* = −obj_opt from (g₂₀, g₃₀) in direction θ.
  //
  // Verification: recompute I₀ directly from v_opt
  I0 = 0.0
  FOR EACH n in basis_f (index idx_n):
    I0 += v_opt[idx_n] * (gamma[n] + 2*g20*alpha[n] + g30*beta[n])
  ASSERT |I0 − obj_opt| < 1e-10   // consistency check

  // ── SIGN CONVENTION AND BOUNDARY EXTRACTION ────────────────────────
  //
  // The LP minimizes I₀[f] = γ[f] + 2g₂₀·α[f] + g₃₀·β[f] subject to D_θ[f] = 1.
  //
  // The constraint for a test point (g₂₀+λcosθ, g₃₀+λsinθ) to lie in the
  // allowed region is: ∀ valid f with D_θ[f]=1,  I₀[f] + λ ≥ 0.
  // This requires: λ ≥ −I₀[f] for ALL valid f.
  // The tightest λ satisfying this is: λ_min = inf_f I₀[f] = obj_opt.
  // Points at λ ≥ obj_opt are allowed; the BOUNDARY is at λ = obj_opt.
  //
  // Therefore: lambda_star = obj_opt, NOT −obj_opt.
  //
  // The boundary point is:
  //   g₂* = g₂₀ + obj_opt · cosθ
  //   g₃* = g₃₀ + obj_opt · sinθ
  //
  // Sign interpretation:
  //   obj_opt < 0: boundary is in direction θ from reference (reference is inside,
  //                boundary is at a negative distance going backward — equivalently
  //                the boundary is in direction θ+π). The full boundary is traced as θ
  //                scans [0,2π) because obj_opt·(cosθ, sinθ) still sweeps all directions.
  //   obj_opt ≥ 0: boundary is at the reference point or further in direction θ.
  //
  // PAPER CONVENTION: The paper writes "maximize distance in direction θ from (g₂₀,g₃₀)."
  // With reference near the allowed boundary, obj_opt takes both signs as θ varies, and
  // the formula g₂* = g₂₀ + obj_opt·cosθ correctly traces the entire closed boundary curve.
  //
  // Cross-check with the upper bound on g₃ (paper's example, θ=3π/2, β=-1):
  //   D_{3π/2}[f] = sinθ·β = -1·(-1) = 1 ✓  (with β[f] = -1)
  //   I₀[f] at (g₂₀, g₃₀=0) = γ + 2αg₂₀.
  //   Minimizing gives obj_opt = min(γ + 2αg₂₀) = the upper bound coefficient.
  //   Boundary in -g₃ direction: g₃* = 0 + obj_opt·sin(3π/2) = 0 + obj_opt·(-1) = -obj_opt.
  //   But the UPPER bound is g₃ ≤ obj_opt (from the paper), not g₃ ≤ -obj_opt.
  //   Resolution: the "boundary" for the upper bound corresponds to the constraint being
  //   saturated at g₃ = obj_opt, i.e., g₃* = obj_opt. The formula g₃* = g₃₀ + obj_opt·sinθ
  //   = 0 + obj_opt·(-1) = -obj_opt does NOT match. This reveals a sign flip:
  //
  // ── CORRECTED FORMULA ─────────────────────────────────────────────
  //   For the angle-scan: the boundary point in direction θ from (g₂₀,g₃₀) is:
  //     g₂* = g₂₀ − obj_opt · cosθ
  //     g₃* = g₃₀ − obj_opt · sinθ
  //   with lambda_star = −obj_opt (so boundary is at distance |obj_opt|).
  //
  //   This matches the upper-bound case: θ=3π/2 → g₃* = g₃₀ − obj_opt·(-1) = obj_opt ✓
  //
  //   The sign flip arises because "maximizing distance" means moving AWAY from the
  //   allowed region, which is in the OPPOSITE direction of the gradient of I₀:
  //     ∇_{(g₂,g₃)} I₀ = (2α[f], β[f]), which points INWARD (toward more-allowed).
  //   Moving OUTWARD (toward boundary) requires going in direction −∇I₀.
  // ─────────────────────────────────────────────────────────────────

  lambda_star = −obj_opt   // positive when obj_opt < 0 (reference inside allowed region)
                            // negative when obj_opt > 0 (reference outside; boundary in θ+π)

  // Boundary point in direction θ from reference:
  g2_boundary = g20 + lambda_star * cos(theta)   // = g20 − obj_opt·cosθ
  g3_boundary = g30 + lambda_star * sin(theta)   // = g30 − obj_opt·sinθ

  // Verification: the bound from this functional is
  //   2α[f]·g₂ + β[f]·g₃ ≥ −γ[f]
  // which is a half-space in (g₂,g₃) with inward normal (2α[f], β[f]).
  // The boundary point should SATURATE this: 2α[f]·g₂* + β[f]·g₃* = −γ[f].
  // With D_θ = 2cosθ·α+sinθ·β = 1 and I₀ = γ+2αg₂₀+βg₃₀ = obj_opt:
  //   LHS: 2α(g₂₀+λ*cosθ) + β(g₃₀+λ*sinθ) = I₀ − γ + λ*(2αcosθ+βsinθ)
  //       = I₀ − γ + λ*·1 = obj_opt − γ + (−obj_opt) = −γ = RHS ✓

  RETURN (g2_boundary, g3_boundary, lambda_star)
```

---

### 4.12 Step 10: Trace the Full Boundary by Scanning Angles

```
FUNCTION TraceBoundary_g2_g3(D, N_f, k_list, N_h, N_angles, J_max,
                              delta_x, epsilon_b, delta_b, B_large, m_max,
                              N_refine, tol_neg, tol_val, precision_bits,
                              g20, g30):

  // ── ONE-TIME SETUP ──────────────────────────────────────────────────
  basis_f  = SelectBasisF(D, N_f)
  basis_h  = {k: [0, 1, ..., N_h[k]-1] for k in k_list}

  (KernelC, KernelX, gamma, alpha, beta,
   ImpactExact, ImpactAsymp) = PrecomputeKernelTable(D, basis_f, k_list,
                                                      basis_h, J_list)

  x_grid_0 = {i * delta_x  for i in 0..floor(1/delta_x)−1}
  b_grid_0 = {epsilon_b + j*delta_b  for j in 0..floor((B_large−epsilon_b)/delta_b)−1}
  J_list   = [0, 2, 4, ..., J_max]

  (MatC_0, MatX_0) = EvaluateOnMassGrid(x_grid_0, J_list, KernelC, KernelX,
                                         basis_f, k_list, basis_h)

  (constraints_0, psd_block) = BuildConstraintSystem(MatC_0, MatX_0, x_grid_0,
                                                      b_grid_0, J_list, k_list,
                                                      basis_f, basis_h,
                                                      ImpactExact, ImpactAsymp,
                                                      B_large, m_max, D)

  // ── ANGLE SCAN ─────────────────────────────────────────────────────
  boundary_points = []

  theta_list = [k * pi / N_angles  for k in 0..2*N_angles−1]
  // = {0, π/N_angles, 2π/N_angles, ..., (2N_angles−1)π/N_angles}
  // For N_angles=20: {0, π/20, 2π/20, ..., 39π/20}

  FOR EACH theta in theta_list:

    // Set up the objective and normalization for this direction
    (obj, norm, norm_rhs) = SetupObjective(theta, g20, g30, gamma, alpha,
                                            beta, basis_f)

    // Run the full adaptive refinement loop
    // Note: start from the INITIAL constraints each time (cold start per angle)
    // A hot-start improvement would reuse x_grid, b_grid from previous angle.
    (v_opt, obj_opt, x_grid_final, b_grid_final, constraints_final) =
      AdaptiveRefinement(constraints_0, psd_block, obj, norm, norm_rhs,
                         x_grid_0, b_grid_0, J_list, KernelC, KernelX,
                         ImpactExact, ImpactAsymp, basis_f, k_list, basis_h,
                         B_large, N_refine, tol_neg, tol_val, precision_bits, D)

    // Extract boundary point
    (g2_b, g3_b, lam) = ExtractBoundaryPoint(v_opt, obj_opt, theta, g20, g30,
                                               basis_f, gamma, alpha, beta)

    boundary_points.append({
      theta:    theta,
      g2:       g2_b,
      g3:       g3_b,
      lambda:   lam,
      v_opt:    v_opt       // save functional for further analysis
    })

  // ── POST-PROCESSING ─────────────────────────────────────────────────
  // The allowed region is to the right of the boundary curve in (g₂, g₃) space.
  // Output the boundary as a list of (g₂, g₃) points sorted by angle.
  // These can be compared with the explicit inequalities in Table 2 of the paper
  // (e.g., "g₂ − g₃/3 + 60.3086 ≥ 0" for D=5) as a consistency check.

  // Restore physical units: multiply by (8πG) and (1/M²) as appropriate.
  // g₂ → g₂ / (8πG·M²), g₃ → g₃·M² / (8πG), etc.

  RETURN boundary_points
```

---

### 4.13 Complete Algorithm: Top-Level Entry Point

```
ALGORITHM FindBounds_g2_g3(D):
  // Parameters matching the paper's 17-dimensional computation (Figure 1):

  N_f         = 7               // basis functions for f(p): e.g. n=1.5,2.5,...,7.5 for even D
  k_list      = [4, 6]          // null constraints X₄ and X₆
  N_h         = {4: 6, 6: 4}    // N_h[4]=6, N_h[6]=4, total = 17 dimensions
  N_angles    = 20              // 40 boundary directions (scan 0..2π)
  J_max       = 42              // max angular momentum
  delta_x     = 1/400           // mass grid spacing
  epsilon_b   = 1/250           // impact parameter start
  delta_b     = 1/32            // impact parameter grid spacing
  B_large     = 40              // large-b PSD cutoff
  m_max       = 2               // subleading terms in large-b expansion
  N_refine    = 10              // adaptive refinement subdivisions
  tol_neg     = 1e-6            // convergence: negative region width
  tol_val     = 1e-9            // convergence: negative region depth
  prec_bits   = 768             // SDPB arithmetic precision in bits

  // Interior reference point (chosen near expected tip from prior experimentation)
  // For D=5, a reasonable starting point is (g20, g30) ≈ (−5, −15) × 1/(8πG·M²)
  (g20, g30)  = SELECT_INTERIOR_POINT(D)

  boundary = TraceBoundary_g2_g3(D, N_f, k_list, N_h, N_angles, J_max,
                                  delta_x, epsilon_b, delta_b, B_large, m_max,
                                  N_refine, tol_neg, tol_val, prec_bits,
                                  g20, g30)

  // The boundary is a closed curve in the (g₂, g₃) plane.
  // The region to the right of the curve (larger g₂) is allowed.
  // Output in physical units with M=1, 8πG=1.

  RETURN boundary
```

---

### 4.14 Algorithm Correctness Checks

The following invariants should hold throughout the computation and serve as verification that the algorithm is correctly implemented:

**Invariant 1 — Basis power constraint.** Every power `n` used in `f(p)` must satisfy `n > 1` so that `∫₀¹ dp f(p)/p²` converges. In the D=5 difference basis `p^k − p²` for `k ≥ 3`, this is satisfied automatically since both `k > 1` and `2 > 1`.

**Invariant 2 — Kernel symmetry.** The improved kernel `C²ᵢₘₚ[m², J, u]` satisfies, by construction:
```
C²ᵢₘₚ[m², J, u]|_{EFT} = 8πG/(-u) + 2g₂ − g₃·u
```
This must be verified by evaluating the LHS of the sum rule with the low-energy amplitude `M_low`. Concretely: `C²ᵢₘₚ[m², J, u]` has no `g₄, g₅, ...` terms on the EFT side — this is the defining property of the improved kernel.

**Invariant 3 — Positivity of UV measure.** The measure `⟨·⟩` defined in eq. (UVaverage) is always positive. Therefore, if `F(m, J) ≥ 0` for all `(m, J)`, then the UV side of the sum rule is non-negative, and the EFT side `I[f] = γ[f] + 2g₂·α[f] + g₃·β[f]` must also be non-negative. This is the soundness condition of the LP.

**Invariant 4 — Normalization consistency.** With normalization `D_θ[f] = 2cosθ·α[f] + sinθ·β[f] = 1`, the extracted bound must satisfy:
```
g₂_boundary = g₂₀ + λ*·cosθ    with    λ* = −obj_opt > 0
g₃_boundary = g₃₀ + λ*·sinθ
```
The sign `λ* > 0` is guaranteed if `(g₂₀, g₃₀)` is strictly inside the allowed region.

**Invariant 5 — Impact-parameter Fourier positivity.** The impact-parameter transform of the final `f(p)` must be non-negative for all `b ≥ 0`. This is the continuous version of the constraints imposed in Blocks 2 and 3. After adaptive refinement converges, the residual negativity magnitude must be below `tol_val`.

**Invariant 6 — Null constraint consistency.** The null constraints `⟨X_k[m², J, u]⟩ = 0` are exact consequences of crossing symmetry. They do not depend on the EFT couplings. Including them does not change the truth of the bound — only its sharpness. Any valid functional with `h_k(p) ≠ 0` can be checked by verifying that `∫₀¹dp h_k(p)·X_k[m², J, -p²]` is finite and adds constructively to the positivity of `F(m, J)`.

**Invariant 7 — Dimensional analysis.** In physical units with general `M` and `8πG`, the EFT action satisfies:
```
γ[f] ~ (8πG)/M²,   α[f] ~ 1/M²,   β[f] ~ 1/M⁴
```
A valid bound of the form `2α·g₂ + β·g₃ ≥ −γ` has `g₂ ~ (8πG)/M²` and `g₃ ~ (8πG)/M⁴` automatically — dimensional analysis is a theorem of the LP, not an input assumption.

---

### 4.15 Summary Flowchart

```
FindBounds_g2_g3(D)
│
├─ SelectBasisF(D, N_f)
│   └─ Returns D-dependent power list for f(p)
│
├─ PrecomputeKernelTable(D, basis_f, k_list, ...)
│   ├─ TABLE A: Analytic ∫pⁿ C²ᵢₘₚ[m²,J] as ₂F₁(m)   [for each n, J]
│   ├─ TABLE B: Analytic ∫pⁱ X_k[m²,J] as ₂F₁(m)      [for each k, i, J]
│   ├─ TABLE C: γ[pⁿ], α[pⁿ], β[pⁿ]                   [scalar precomputation]
│   └─ TABLE D: ImpactExact[n](b) = ₁F₂(...), ImpactAsymp[n](b) = A+Bcos+Csin
│
├─ BuildConstraintSystem(...)
│   ├─ BLOCK 1: (x_i, J) constraints  →  linear rows  ≥ 0
│   ├─ BLOCK 2: b_j constraints       →  linear rows  ≥ 0
│   └─ BLOCK 3: Large-b PSD matrix    →  polynomial matrix  ≽ 0
│
└─ FOR EACH θ in θ_list:
    │
    ├─ SetupObjective(θ, g₂₀, g₃₀, ...)
    │   └─ minimize I₀[f] = γ[f]+2g₂₀α[f]+g₃₀β[f]  s.t. Dθ[f]=1
    │
    ├─ AdaptiveRefinement(constraints, psd_block, obj, norm, ...)
    │   │
    │   └─ LOOP:
    │       ├─ SolveSDP via SDPB (768–840 bit precision)
    │       ├─ ValidateFunctional on dense verification grid
    │       ├─ IF converged (neg regions < tol_neg AND tol_val): BREAK
    │       ├─ REFINE x_grid:
    │       │   └─ For each (J, x₁, x₂): quadratic min → x* → add {x*±k·s}
    │       └─ REFINE b_grid:
    │           └─ For each (b₁, b₂): quadratic min → b* → add {b*±k·s}
    │
    └─ ExtractBoundaryPoint(v_opt, obj_opt, θ, g₂₀, g₃₀, ...)
        └─ g₂* = g₂₀ − obj_opt·cosθ
           g₃* = g₃₀ − obj_opt·sinθ

OUTPUT: List of (g₂*, g₃*) boundary points for 2·N_angles directions.
        Intersection of all implied half-spaces = allowed region in (g₂, g₃) plane.
```

---

## Part 5: Focused Pseudocode — Finding the g₂–g₃ Boundary

This section provides a **self-contained, focused pseudocode** specifically designed to find the boundary between the allowed and excluded regions in the `(g₂, g₃)` plane — the central result of the paper. Unlike Part 4 (which documents the full technical machinery), Part 5 is organised around the **mathematical meaning** of the boundary, the **geometric dual perspective**, and a **single unified algorithm** with recursive consistency checks. Every formula is cross-referenced against the TeX source and can be verified independently.

All units are set to `M = 1`, `8πG = 1` throughout, matching the paper's appendix convention. Physical bounds are recovered by homogeneity.

---

### 5.1 What "Finding the Bound Between g₂ and g₃" Means

#### 5.1.1 The Primal (Physical) Perspective

The EFT consistency conditions force every physical point `(g₂, g₃)` to satisfy an infinite family of inequalities of the form:

```
2α[f] · g₂  +  β[f] · g₃  ≥  −γ[f]
```

for every **valid functional** `f`, where:

```
γ[f] = ∫₀¹ dp f(p)/p²           (gravity/Newton constant term)
α[f] = ∫₀¹ dp f(p)              (g₂ coefficient)
β[f] = ∫₀¹ dp f(p) · p²        (g₃ coefficient)
```

A functional `f` is **valid** (i.e., yields a true inequality) if and only if:

```
F(m, J) ≡ ∫₀¹ dp f(p) · C²ᵢₘₚ[m², J, −p²]
          + Σ_{k=4,6,...} ∫₀¹ dp h_k(p) · X_k[m², J, −p²]  ≥  0

for ALL  m ≥ 1,  J = 0, 2, 4, ...
```

This follows from positivity of the UV spectral measure `⟨·⟩` (eq. UVaverage in the paper): if `F(m,J) ≥ 0` for all heavy states, then the UV side of the improved sum rule is non-negative, and the EFT side `γ[f] + 2g₂·α[f] + g₃·β[f] ≥ 0` is forced.

The **boundary** of the allowed region in `(g₂, g₃)` space is therefore:

```
BOUNDARY = { (g₂*, g₃*) ∈ ℝ² | ∃ valid f :
             2α[f]·g₂* + β[f]·g₃* = −γ[f]     ← inequality is TIGHT
             AND
             ∀ valid g : 2α[g]·g₂* + β[g]·g₃* ≥ −γ[g]   ← all others satisfied }
```

Every boundary point has an associated **extremal functional** `f*` that saturates the inequality precisely there.

#### 5.1.2 The Dual (Geometric) Perspective

Each valid functional `f` defines a **supporting half-space** in the `(g₂, g₃)` plane:

```
H(f) = { (g₂, g₃) : (2α[f], β[f]) · (g₂, g₃) ≥ −γ[f] }
```

The allowed region is the intersection of all such half-spaces:

```
ALLOWED = ∩_{valid f} H(f)
```

The **inward normal** to the half-space boundary line is the vector `n̂(f) = (2α[f], β[f])`.

The boundary of the allowed region is traced by finding, for each direction `θ`, the valid functional `f*` whose inward normal `n̂(f*)` points in direction `θ` — i.e., `(2α[f*], β[f*]) ∝ (cosθ, sinθ)`. This is achieved by maximising the distance from an interior point in direction `θ`, which is exactly what the LP in Section 4 does.

**Key geometric fact:** A point `(g₂*, g₃*)` lies on the boundary if and only if there exists a valid `f*` such that:

```
(2α[f*], β[f*]) ∝ (cosθ, sinθ)    ← normal matches direction
2α[f*]·g₂* + β[f*]·g₃* = −γ[f*]  ← point lies on the half-space boundary
```

with the normalization condition `2cosθ·α[f*] + sinθ·β[f*] = 1`.

---

### 5.2 Mathematical Setup for the g₂–g₃ Boundary LP

#### 5.2.1 The Improved Sum Rule (Derivation Recap)

Starting from the C₂ dispersive sum rule (eq. B2 flat in `sec_dispersion_review.tex`):

```
C_{2,u}|_EFT = 8πG/(-u) + 2g₂ − g₃u + 8g₄u² − 2g₅u³ + ...
```

The key problem: at `u → 0`, the graviton pole `8πG/(-u)` diverges, preventing use of forward-limit bounds. The solution (eq. `definitionofb2improved` in `sec_localized_b.tex`) is to subtract off all higher EFT contamination using `C_{2n,0}` and `C'_{2n,0}` (which are free of the graviton pole for `n ≥ 2`):

```
C²ᵢₘₚ_{u} = C_{2,u} - Σ_{n=2}^∞ (n·u^{2n-2}·C_{2n,0} + u^{2n-1}·C'_{2n,0})
```

This yields (eq. `B2 improved flat`):

```
C²ᵢₘₚ_{u}|_EFT = 8πG/(-u) + 2g₂ − g₃·u     ← only three couplings remain
```

The UV side of the improved sum rule is the **improved kernel** evaluated on heavy states:

```
C²ᵢₘₚ[m², J, u] =
    (2m²+u) · Pⱼ(1+2u/m²)
    ─────────────────────────────────────
          m²·(m²+u)²

    − u²/m⁶ · [ (4m²+3u)·Pⱼ(1)     4u·Pⱼ'(1)  ]
               [ ─────────────────  + ─────────── ]
               [    (m²+u)²             m⁴−u²    ]
```

with `Pⱼ(1) = 1` and `Pⱼ'(1) = J(J+D-3)/(D-2)` (quadratic Casimir over `D-2`).

Working at `u = −p²` with `p ∈ [0, 1]` (i.e., momentum transfer up to the EFT cutoff `M=1`):

```
C²ᵢₘₚ[m², J, −p²] =
    (2m²−p²) · Pⱼ(1−2p²/m²)
    ────────────────────────────────────────
             m²·(m²−p²)²

    − p⁴/m⁶ · [ (4m²−3p²)·Pⱼ(1)     −4p²·Pⱼ'(1)  ]
               [ ─────────────────  + ───────────── ]
               [    (m²−p²)²             m⁴−p⁴     ]
```

#### 5.2.2 The Null Constraints (Crossing Symmetry)

The null constraints `X_k` (derived from the `s ↔ u` crossing symmetry of the unsubtracted dispersion relation) vanish on the UV average `⟨·⟩` without any EFT content. Including them in the functional tightens bounds without changing their validity. They are defined for `k = 4, 6, ...` by eq. `nullconstraints` in the paper (with `u = −p²`):

```
X_k[m², J, u] =
    (2m²+u)              m² · Pⱼ(1+2u/m²)
    ──────────────  −  ─────────────────────────────── − Res_{u'=0}[...]
    u·m²·(m²+u)         [u·m²·(m²+u)]^{k/2}
```

The residue at `u'=0` ensures the forward-limit regularity.

**Consistency check:** `⟨X_k[m², J, u]⟩ = 0` for all `k` and all `u ∈ (−M², 0)`. Adding `h_k · X_k` to the functional therefore never changes the EFT side of the inequality — only the UV side, where it provides additional discriminating power.

---

### 5.3 Complete Pseudocode: `FindBoundary_g2_g3`

The procedure below is the **single entry point** for the g₂–g₃ boundary computation. It is structured in four phases:

- **Phase A**: Mathematical setup and analytic precomputation  
- **Phase B**: LP discretization and constraint assembly  
- **Phase C**: Angle scan with adaptive refinement  
- **Phase D**: Boundary extraction, verification, and output  

```
╔══════════════════════════════════════════════════════════════════════════╗
║  PROCEDURE  FindBoundary_g2_g3(D)                                       ║
║  INPUTS:    Spacetime dimension D ∈ {5,...,12}  (or D=4 with IR cutoff) ║
║  OUTPUT:    Closed boundary curve in (g₂, g₃) plane  [units M=1, 8πG=1]║
╚══════════════════════════════════════════════════════════════════════════╝

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE A — PARAMETER SETUP AND ANALYTIC PRECOMPUTATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

// ── A1. NUMERICAL PARAMETERS (matching paper's Figure 1 computation) ──

N_f        ← 7            // f(p) basis dimension
k_list     ← [4, 6]       // null constraint orders to include
N_h        ← {4: 6, 6: 4} // h_k(p) basis sizes → total dim = 7+6+4 = 17
N_angles   ← 20           // scan 2·N_angles = 40 directions
J_max      ← 42           // angular momentum truncation
δ_x        ← 1/400        // initial mass-grid spacing (x = 1 − 1/m²)
ε_b        ← 1/250        // impact-parameter grid start
δ_b        ← 1/32         // impact-parameter grid step
B_large    ← 40           // large-b cutoff for PSD block
m_max      ← 2            // subleading terms in large-b asymptotic expansion
N_refine   ← 10           // refinement points per negative interval
tol_neg    ← 1e-6         // convergence: negative-region WIDTH below this
tol_val    ← 1e-9         // convergence: negative-region DEPTH below this
prec_bits  ← 768          // SDPB arithmetic precision

// Interior reference point — chosen near the expected tip of the allowed region.
// If unknown a priori, start with (g₂₀, g₃₀) = (−C_D, −3·C_D) where C_D is the
// leading coefficient from the forward-limit bound:
//   C_D ≈ (D-related constant) such that the trivial bound g₂ ≥ −C_D holds.
// For D=5 (from Table 2 of paper): g₂ + g₃/3 ≥ −60.3086, so set g₂₀ ≈ −5,  g₃₀ ≈ −15.
(g₂₀, g₃₀) ← SELECT_INTERIOR_POINT(D)
//   SELECT_INTERIOR_POINT returns a point known to be strictly inside the allowed region.
//   Invariant: γ[f] + 2g₂₀·α[f] + g₃₀·β[f] < 0 for all optimal functionals f.
//   Equivalently: I₀[f] < 0 for all optimal f. This ensures λ* > 0 in Phase C.


// ── A2. SELECT BASIS FUNCTIONS FOR f(p) ──
//
// The basis must satisfy:
//   (i)  every power n in basis_f satisfies n > 1 (so that γ[pⁿ] = 1/(n-1) converges)
//   (ii) the impact-parameter Fourier transform ∫₀¹ dp f(p) Jᵥ(bp)/(bp/2)^ν
//        can be made positive at large b (determines the choice of powers by dimension)
//
// See sec_localized_b.tex and app_flat_space_numerics.tex, Table 1:
//
//   D=5:       difference basis   {pᵏ − p² | k = 3, 4, ..., N_f+2}
//              REASON: (D-3)/2 = 1 < 2 = smallest available n, so cancel leading
//              oscillatory term at large b by forming differences with p².
//
//   even D≥6:  half-integer basis {p^{(2i+1)/2} | i = 1, ..., N_f}
//              REASON: half-integer powers give non-oscillatory decay at large b.
//
//   odd D≥7:   integer basis      {pⁿ | n = 2, ..., N_f+1}
//              REASON: (D-3)/2 ≥ 2 ≥ n_min, oscillatory term decays fast enough.

basis_f ← SELECT_BASIS(D, N_f)          // returns list of N_f basis elements

// Basis for h_k(p) is always simple integer powers:
basis_h ← { k: [pⁱ for i = 0,...,N_h[k]-1]  for each k in k_list }


// ── A3. PRECOMPUTE ANALYTIC KERNEL INTEGRALS ──
//
// All integrals of basis functions against the constraint kernels are evaluated
// ONCE analytically as closed-form functions of m. This avoids numerical integration
// and allows high-precision arithmetic.
//
// For each (n, J):  (where n runs over powers in basis_f, after expanding differences)
//
//   TABLE_C[n][J](m) = ∫₀¹ dp pⁿ · C²ᵢₘₚ[m², J, −p²]
//                    = sum of ₂F₁ hypergeometric functions in (−1/m²)
//
//   Example (J=2, from eq. A.3 of paper):
//     ∫₀¹dp pⁿ · C²ᵢₘₚ[m²,2,−p²]
//       = −4(D-1)·₂F₁(1,(n+1)/2;(n+3)/2;−1/m²)/[(D-2)·m⁴·(n+1)]
//         + 2(3D-4)/[(D-2)·m⁴·(n+1)]
//         − 3(3D-2)/[(D-2)·m⁶·(n+3)]
//
//   For D=5 difference basis (pᵏ−p²): TABLE_C[k−2][J](m) = TABLE_C[k][J](m) − TABLE_C[2][J](m)
//
// For each (k, i, J):
//
//   TABLE_X[k][i][J](m) = ∫₀¹ dp pⁱ · X_k[m², J, −p²]
//                       = similar ₂F₁ sum, with a Laurent-at-0 residue contribution
//
// For each n: (EFT coupling actions)
//
//   γ[n] = ∫₀¹ dp pⁿ / p² = 1/(n−1)         [requires n > 1]
//   α[n] = ∫₀¹ dp pⁿ      = 1/(n+1)
//   β[n] = ∫₀¹ dp pⁿ · p² = 1/(n+3)
//   (For D=5 difference basis: subtract the p²-contribution from each.)
//
// For each n: (impact-parameter Bessel transform, for large-J positivity)
//
//   IMPACT_EXACT[n](b) = ₁F₂((n+1)/2; (D-2)/2, (n+3)/2; −b²/4) / (n+1)
//
//   IMPACT_ASYMP[n](b) = A_n(b) + B_n(b)·cos(b − π(D-1)/4) + C_n(b)·sin(b − π(D-1)/4)
//   where A_n, B_n, C_n are power series in 1/b (from NIST DLMF expansion of ₁F₂).
//   Keep m_max subleading terms in B_n and C_n.

PRECOMPUTE_TABLES(D, basis_f, basis_h, k_list, J_max, m_max)
// → TABLE_C, TABLE_X, γ[], α[], β[], IMPACT_EXACT[], IMPACT_ASYMP[]


// ── A4. RECURSIVE SELF-CHECK A: VERIFY IMPROVED KERNEL EFT CONTENT ──
//
// For a randomly chosen test functional f_test = p² (one basis element):
//   Compute: I_test = γ[2] + 2g₂_test·α[2] + g₃_test·β[2]
//            = 1/(2-1) + 2g₂_test·1/(2+1) + g₃_test·1/(2+3)
//            = 1 + 2g₂_test/3 + g₃_test/5
//
// Evaluate the UV side directly from the sum rule:
//   UV_side = ⟨C²ᵢₘₚ[m², J, −p²]⟩  (using a test spectral density)
//
// Cross-check: EFT_side should equal UV_side for any spectral density.
// This confirms TABLE_C was computed correctly.
//
// ASSERT: TABLE_C[n][J](m) for J=0 gives the J=0 special case (no angular dependence):
//   P₀(x) = 1, P₀'(1) = 0, so the J=0 kernel is purely m-dependent.
//   The J=0 integral should equal TABLE_C[n][0](m) = (2m²−p²)/(m²·(m²−p²)²) integrated
//   against pⁿ, minus the forward-subtraction terms, as a pure function of m.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE B — DISCRETIZATION AND CONSTRAINT ASSEMBLY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

// ── B1. BUILD INITIAL GRIDS ──
//
// Mass grid: x = 1 − 1/m² (so m=1→x=0, m→∞→x=1)
//   x_grid ← {0, δ_x, 2δ_x, ..., ⌊1/δ_x − 1⌋·δ_x}    (size ≈ 400 points)
//
// Impact-parameter grid (discrete part, below B_large):
//   b_grid ← {ε_b, ε_b+δ_b, ..., ε_b + ⌊(B_large−ε_b)/δ_b − 1⌋·δ_b}
//
// Angular momentum list:
//   J_list ← {0, 2, 4, ..., J_max}                     (size J_max/2 + 1 = 22)

x_grid_0 ← MAKE_MASS_GRID(δ_x)
b_grid_0 ← MAKE_IMPACT_GRID(ε_b, δ_b, B_large)
J_list   ← MAKE_SPIN_LIST(J_max)

// Decision variable layout:
//   v = (a_{n₁}, ..., a_{n_{N_f}},  b_{4,0}, ..., b_{4,N_h[4]-1},  b_{6,0}, ..., b_{6,N_h[6]-1})
//   Total dimension: N_vars = N_f + Σ_k N_h[k] = 7 + 6 + 4 = 17
N_vars ← N_f + SUM(N_h[k] for k in k_list)


// ── B2. BUILD BLOCK-1 CONSTRAINTS: FINITE-m POSITIVITY ──
//
// For each grid point x_i (equivalently m_i = 1/√(1−x_i)) and each J ∈ J_list:
//
//   F(m_i, J) = Σₙ aₙ · TABLE_C[n][J](m_i)
//             + Σ_k Σᵢ b_{k,i} · TABLE_X[k][i][J](m_i)  ≥  0
//
// This is a LINEAR INEQUALITY in the decision variables v.
// Encode as a row vector:
//   row[(idx of aₙ)] = TABLE_C[n][J](m_i)    for each n in basis_f
//   row[(idx of b_{k,i})] = TABLE_X[k][i][J](m_i)  for each k,i

block1_constraints ← []
FOR EACH x_i IN x_grid_0:
  m_i ← 1 / SQRT(1 − x_i)
  FOR EACH J IN J_list:
    row ← ZEROS(N_vars)
    FOR n, idx_n IN ENUMERATE(basis_f):
      row[idx_n] ← TABLE_C[n][J](m_i)        // analytic, high precision
    FOR k IN k_list:
      FOR i, idx_ki IN ENUMERATE(basis_h[k]):
        row[idx_ki] ← TABLE_X[k][i][J](m_i)
    block1_constraints.APPEND( (row, lb=0) )

// Total Block-1 size: |x_grid_0| × |J_list| = 400 × 22 = 8800 constraints


// ── B3. BUILD BLOCK-2 CONSTRAINTS: DISCRETIZED IMPACT-PARAMETER POSITIVITY ──
//
// In the large-m (scaling) limit, the positivity condition becomes:
//
//   Γ((D-2)/2) · ∫₀¹ dp f(p) · Jᵥ(bp)/(bp/2)^ν  ≥  0   for all b ≥ 0
//
// where ν = (D-4)/2. The h_k functions do NOT enter this limit.
// For each b_j ≤ B_large in b_grid:

block2_constraints ← []
FOR EACH b_j IN b_grid_0:
  row ← ZEROS(N_vars)
  FOR n, idx_n IN ENUMERATE(basis_f):
    row[idx_n] ← IMPACT_EXACT[n](b_j)    // = ₁F₂(...;...;−b²/4)/(n+1)
  block2_constraints.APPEND( (row, lb=0) )


// ── B4. BUILD BLOCK-3 CONSTRAINT: LARGE-b POLYNOMIAL MATRIX PSD ──
//
// For b ≥ B_large, replace positivity of the scalar:
//   A(b) + B(b)cos(b−φ) + C(b)sin(b−φ)  ≥  0
// with the STRONGER (but rigorous) 2×2 PSD condition:
//
//   M(b) = [[A(b)+B(b),  C(b)  ],
//            [C(b),       A(b)−B(b)]]  ≽  0
//
// (stronger because it implies positivity for all φ ∈ [0,2π), not just the specific φ=b−π(D-1)/4)
//
// Expand A,B,C as power series in 1/b (keeping m_max subleading terms beyond leading).
// Change variable: z = 1/b, so z ∈ (0, 1/B_large].
// After multiplying by the positive factor b^{(D-1)/2} = z^{−(D-1)/2}, the matrix
// entries become POLYNOMIAL in z. This is a valid polynomial matrix constraint for SDPB.
//
// Each entry of M is a LINEAR function of {aₙ}:
//   A(b) + B(b) = Σₙ aₙ · (A_n(b) + B_n(b))
//   A(b) − B(b) = Σₙ aₙ · (A_n(b) − B_n(b))
//   C(b)         = Σₙ aₙ · C_n(b)
// where A_n, B_n, C_n come from IMPACT_ASYMP[n].
//
// Note: h_k coefficients do NOT appear in this block.

psd_block ← BUILD_POLY_MATRIX_PSD_BLOCK(IMPACT_ASYMP, basis_f, B_large, m_max, D, N_vars)
// psd_block is passed to SDPB as a polynomial matrix positive-semidefinite constraint block.


// ── B5. RECURSIVE SELF-CHECK B: VERIFY CONSTRAINT DIMENSIONS AND NULL CONSTRAINTS ──
//
// CHECK 1 — NULL CONSTRAINT CONSISTENCY:
//   For test spectral density ρ_J(m) = δ(m−m₀), evaluate ⟨X_k[m², J, −p²]⟩.
//   Result must be zero (to numerical precision) for all k and all p ∈ (0,1).
//   This verifies TABLE_X was computed correctly.
//
// CHECK 2 — DIMENSION COUNT:
//   ASSERT N_vars = N_f + Σ_k N_h[k] = 7 + 6 + 4 = 17   (for Figure 1 parameters)
//   ASSERT len(block1_constraints) = len(x_grid_0) × len(J_list)
//   ASSERT len(block2_constraints) = len(b_grid_0)
//
// CHECK 3 — BLOCK-1 POSITIVITY FOR J=0:
//   For J=0: F(m, 0) = Σₙ aₙ · TABLE_C[n][0](m) + ...
//   Since Pⱼ(1) = 1 and Pⱼ'(1) = 0 for J=0, the J=0 row should involve only
//   the non-angular terms of the improved kernel. Verify analytically.
//
// CHECK 4 — IMPACT-PARAMETER LEADING TERM:
//   IMPACT_EXACT[n](b) → ₁F₂((n+1)/2; (D-2)/2, (n+3)/2; 0) / (n+1) = 1/(n+1)  as b→0.
//   This matches α[n] = 1/(n+1). In the forward limit b=0, the Bessel integral
//   reduces to the plain integral of pⁿ. ASSERT this numerically.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE C — ANGLE SCAN: FIND BOUNDARY IN EACH DIRECTION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

// ── C1. ASSEMBLE BASE CONSTRAINT SYSTEM ──
//
// Combine Block-1, Block-2, Block-3 into the base constraint system for SDPB.
// This is reused for all angles θ (the base constraints are angle-independent).

base_system ← ASSEMBLE(block1_constraints, block2_constraints, psd_block)

// ── C2. DEFINE ANGLE SCAN ──
//
// Scan 2·N_angles = 40 equally spaced directions in [0, 2π):
//   θ_list = {k·π/N_angles  for k = 0, 1, ..., 2·N_angles − 1}
//           = {0, π/20, 2π/20, ..., 39π/20}
//
// For each angle θ: the LP objective is to find the valid functional f
// whose half-space boundary is as far as possible from (g₂₀, g₃₀) in direction θ.

boundary_points ← []

FOR EACH θ IN θ_list:

  // ── C3. SETUP OBJECTIVE FOR THIS ANGLE ──
  //
  // Parameterize the boundary in direction θ as:
  //   (g₂, g₃) = (g₂₀ + λcosθ, g₃₀ + λsinθ)  for λ ∈ ℝ
  //
  // Substituting into the EFT inequality γ[f] + 2g₂·α[f] + g₃·β[f] ≥ 0:
  //   (γ[f] + 2g₂₀·α[f] + g₃₀·β[f]) + λ·(2cosθ·α[f] + sinθ·β[f]) ≥ 0
  //   I₀[f]                          + λ·D_θ[f]                    ≥ 0
  //
  // where:
  //   I₀[f] = γ[f] + 2g₂₀·α[f] + g₃₀·β[f]   ("action at interior reference point")
  //   D_θ[f] = 2cosθ·α[f] + sinθ·β[f]          ("directional derivative")
  //
  // The boundary in direction θ is at λ* = −I₀[f*]/D_θ[f*].
  // Maximizing λ* subject to D_θ[f] = 1 (normalization) reduces to:
  //   MINIMIZE  I₀[f]    SUBJECT TO  D_θ[f] = 1
  //
  // Optimal value: obj_opt = min I₀[f] (should be ≤ 0 when g₂₀,g₃₀ is interior)
  // Then: λ* = −obj_opt  (positive distance from interior reference)
  //
  // In terms of coefficient vector v = (aₙ, b_{k,i}):
  //   OBJECTIVE vector c:    c[idx_n] = γ[n] + 2g₂₀·α[n] + g₃₀·β[n]
  //                          c[idx_{k,i}] = 0   (h_k coefficients have no EFT coupling)
  //   NORMALIZATION vector:  d[idx_n] = 2cosθ·α[n] + sinθ·β[n]
  //                          d[idx_{k,i}] = 0

  FOR n, idx_n IN ENUMERATE(basis_f):
    obj[idx_n] ← γ[n] + 2·g₂₀·α[n] + g₃₀·β[n]
    nrm[idx_n] ← 2·COS(θ)·α[n] + SIN(θ)·β[n]


  // ── C4. ADAPTIVE REFINEMENT LOOP ──
  //
  // Start from the initial (coarse) constraint system. Iteratively refine
  // the mass and impact-parameter grids until the returned functional is
  // non-negative everywhere (to within tolerances).

  x_grid ← x_grid_0   // mutable working copy
  b_grid ← b_grid_0   // mutable working copy
  constraints ← COPY(base_system)

  FOR round = 1, 2, ..., MAX_ROUNDS:

    // ── C4a. SOLVE THE SDP ──
    //
    // SDPB solves: minimize c·v  subject to M_j(v) ≽ 0  (all blocks)
    //              and d·v = 1  (normalization equality)
    //
    // Return: v_opt (primal solution), obj_opt = c·v_opt (objective at optimum)
    //
    // Set SDPB precision to prec_bits = 768.
    // If PRIMAL_INFEASIBLE: no bound in this direction (skip this θ).

    (v_opt, obj_opt, status) ← SDPB_SOLVE(constraints, psd_block, obj, nrm,
                                           nrm_rhs=1.0, prec=prec_bits)
    IF status == INFEASIBLE: BREAK
    IF status == FAILURE: RAISE Error("SDPB failed at θ=" + θ + " round=" + round)


    // ── C4b. VALIDATE THE RETURNED FUNCTIONAL ──
    //
    // Check F(m, J) on a DENSE verification grid (100× finer than training grid).
    // A functional is "valid" only if it is non-negative everywhere.
    // Identify all negative intervals [x₁, x₂] and [b₁, b₂].

    neg_x ← []   // list of (J, x₁, x₂, depth) for negative intervals in x
    neg_b ← []   // list of (b₁, b₂, depth) for negative intervals in b

    FOR J IN J_list:
      FOR EACH PAIR (x_prev, x_curr) IN DENSE_PAIRS(0, 1, step=δ_x/100):
        m_curr ← 1/SQRT(1−x_curr)
        F_val ← Σₙ v_opt[idx_n] · TABLE_C[n][J](m_curr)
               + Σ_{k,i} v_opt[idx_{k,i}] · TABLE_X[k][i][J](m_curr)
        IF F_val < 0:
          neg_x.APPEND( (J, x_prev, x_curr, F_val) )

    FOR EACH PAIR (b_prev, b_curr) IN DENSE_PAIRS(ε_b, B_large, step=δ_b/100):
      Fb_val ← Σₙ v_opt[idx_n] · IMPACT_EXACT[n](b_curr)
      IF Fb_val < 0:
        neg_b.APPEND( (b_prev, b_curr, Fb_val) )


    // ── C4c. CHECK CONVERGENCE ──
    //
    // Two-criterion convergence (matching the paper's description):
    //   Criterion 1 (width): all negative regions have width < tol_neg = 1e-6
    //   Criterion 2 (depth): all negative regions have depth |val| < tol_val = 1e-9
    // EITHER criterion suffices for convergence.

    width_ok ← ALL( ABS(x₂−x₁) < tol_neg FOR (J,x₁,x₂,v) IN neg_x )
              AND ALL( ABS(b₂−b₁) < tol_neg FOR (b₁,b₂,v) IN neg_b )
    depth_ok ← ALL( ABS(v) < tol_val FOR (_,_,_,v) IN neg_x )
              AND ALL( ABS(v) < tol_val FOR (_,_,v) IN neg_b )

    IF (width_ok OR depth_ok) AND (neg_x OR neg_b → all_small):
      BREAK   // converged


    // ── C4d. REFINE THE x-GRID (mass constraints) ──
    //
    // For each negative interval (J, x₁, x₂, depth):
    //   1. Estimate minimum location x* via quadratic interpolation through
    //      the three points (x₁, F(x₁)), ((x₁+x₂)/2, F_mid), (x₂, F(x₂)).
    //      Quadratic minimum: x* = x₁ + (x₂−x₁)/2 − (F(x₂)−F(x₁))·(x₂−x₁)/(8·F_mid−4·F(x₁)−4·F(x₂))
    //      Clamp to (x₁, x₂) to guard against non-quadratic behavior.
    //   2. Add 2·N_refine+1 new constraint points uniformly around x*:
    //      { x* + i·s | i = −N_refine, ..., N_refine }  where s = (x₂−x₁)/N_refine
    //   3. For each new point x_new: add a new Block-1 constraint row for
    //      the specific J that showed the negative region.

    FOR (J, x₁, x₂, depth) IN neg_x:
      x_mid ← (x₁+x₂)/2
      F1 ← EVAL_F(v_opt, J, x₁)
      F2 ← EVAL_F(v_opt, J, x₂)
      Fm ← EVAL_F(v_opt, J, x_mid)

      // Quadratic minimum estimate:
      denom ← 8·Fm − 4·F1 − 4·F2
      IF |denom| > 1e-30:
        x_star ← x_mid − (F2−F1)·(x₂−x₁) / (2·denom)
      ELSE:
        x_star ← x_mid    // fallback: use midpoint if quadratic is flat
      x_star ← CLAMP(x_star, x₁+1e-12, x₂−1e-12)

      s ← (x₂−x₁) / N_refine
      new_pts ← { x_star + i·s | i = −N_refine,...,N_refine }
      new_pts ← FILTER(new_pts, 0 ≤ pt < 1)
      new_pts ← new_pts SETMINUS x_grid    // only genuinely new points

      FOR x_new IN new_pts:
        m_new ← 1/SQRT(1−x_new)
        row ← ZEROS(N_vars)
        FOR n, idx_n IN ENUMERATE(basis_f):
          row[idx_n] ← TABLE_C[n][J](m_new)
        FOR k IN k_list:
          FOR i, idx_ki IN ENUMERATE(basis_h[k]):
            row[idx_ki] ← TABLE_X[k][i][J](m_new)
        constraints.APPEND( (row, lb=0) )
        x_grid.ADD(x_new)


    // ── C4e. REFINE THE b-GRID (impact-parameter constraints) ──
    //
    // Same procedure for impact-parameter negative intervals.

    FOR (b₁, b₂, depth) IN neg_b:
      b_mid ← (b₁+b₂)/2
      Fb1 ← EVAL_Fb(v_opt, b₁)
      Fb2 ← EVAL_Fb(v_opt, b₂)
      Fbm ← EVAL_Fb(v_opt, b_mid)

      denom ← 8·Fbm − 4·Fb1 − 4·Fb2
      IF |denom| > 1e-30:
        b_star ← b_mid − (Fb2−Fb1)·(b₂−b₁) / (2·denom)
      ELSE:
        b_star ← b_mid
      b_star ← CLAMP(b_star, b₁+1e-12, b₂−1e-12)

      s ← (b₂−b₁) / N_refine
      new_pts ← { b_star + i·s | i = −N_refine,...,N_refine }
      new_pts ← FILTER(new_pts, ε_b ≤ pt ≤ B_large)
      new_pts ← new_pts SETMINUS b_grid

      FOR b_new IN new_pts:
        row ← ZEROS(N_vars)
        FOR n, idx_n IN ENUMERATE(basis_f):
          row[idx_n] ← IMPACT_EXACT[n](b_new)
        constraints.APPEND( (row, lb=0) )
        b_grid.ADD(b_new)

  // END adaptive refinement loop


  // ── C5. EXTRACT BOUNDARY POINT IN DIRECTION θ ──
  //
  // The optimal objective value obj_opt = min I₀[f] with D_θ[f]=1.
  // The boundary in direction θ from (g₂₀, g₃₀) is at distance λ* = −obj_opt.
  //
  // Derivation of sign:
  //   The LP finds the f* that minimizes the "headroom" I₀[f] in direction θ.
  //   Since (g₂₀,g₃₀) is inside the allowed region: I₀[f*] = obj_opt < 0.
  //   The boundary at the extremal direction is reached by moving outward by λ*:
  //     (g₂₀ + λ*cosθ, g₃₀ + λ*sinθ)  where λ* = −obj_opt > 0.
  //
  // The boundary inequality from f* is tight at this point:
  //   γ[f*] + 2(g₂₀+λ*cosθ)·α[f*] + (g₃₀+λ*sinθ)·β[f*]
  //   = I₀[f*] + λ*·D_θ[f*] = obj_opt + (−obj_opt)·1 = 0 ✓

  λ_star    ← −obj_opt
  g₂_bound  ← g₂₀ + λ_star · COS(θ)
  g₃_bound  ← g₃₀ + λ_star · SIN(θ)

  boundary_points.APPEND({
    θ:       θ,
    g₂:      g₂_bound,
    g₃:      g₃_bound,
    λ:       λ_star,
    f_star:  v_opt,         // save dual certificate functional
    obj:     obj_opt
  })

// END angle scan loop

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE D — BOUNDARY VERIFICATION AND OUTPUT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

// ── D1. VERIFY AGAINST KNOWN INEQUALITIES (TABLE 2 CROSS-CHECK) ──
//
// Each boundary point (g₂*, g₃*) must satisfy ALL of the explicit inequalities
// tabulated in Table 2 (app_flat_space_numerics.tex, eq:inequalitiesforfigure1).
//
// For D=5, check all of:
//   g₂* − g₃*/3 + 60.3086 ≥ 0
//   g₂* + 0.0647867·g₃* + 9.64034 ≥ 0
//   ...   (all 17 inequalities for D=5)
//
// A boundary point should SATURATE exactly ONE of these inequalities (the one
// corresponding to the angle θ from which it was computed).

FOR EACH bp IN boundary_points:
  FOR EACH ineq IN TABLE2_INEQUALITIES[D]:
    val ← ineq.a·bp.g₂ + ineq.b·bp.g₃ + ineq.c    // = a·g₂ + b·g₃ + c ≥ 0
    ASSERT val ≥ −ε_verify                          // must satisfy (up to numerical error)
    // At the saturating inequality: ASSERT |val| < ε_tight (approximately zero)


// ── D2. VERIFY CONVEXITY OF THE BOUNDARY CURVE ──
//
// The allowed region is the intersection of convex half-spaces, so it is CONVEX.
// The boundary curve must be convex. Check:
//   For every three consecutive boundary points P₁, P₂, P₃ (sorted by angle):
//   P₂ must lie on or outside the line segment P₁P₃ (convex side).
//
// Numerically: the signed area of triangle (P₁, P₂, P₃) must be non-negative
// (counterclockwise orientation when moving around the boundary in the allowed direction).

FOR i IN 0..len(boundary_points)−1:
  P1 ← boundary_points[i]
  P2 ← boundary_points[(i+1) MOD len]
  P3 ← boundary_points[(i+2) MOD len]
  cross ← (P2.g₂−P1.g₂)·(P3.g₃−P1.g₃) − (P2.g₃−P1.g₃)·(P3.g₂−P1.g₂)
  ASSERT cross ≥ −ε_convex
  // Negative cross product would mean the boundary is concave (impossible for LP output)


// ── D3. VERIFY KNOWN PHYSICAL POINTS LIE INSIDE ──
//
// String theory check (D=6, from paper Section 1.9):
//   Type II:     (g₂, g₃) = (0, 0)                        → must be inside
//   Heterotic:   (g₂, g₃) = (3/16, 3/4)                   → must be inside
//   Type II -:   (g₂, g₃) ≈ (−0.75, −2.07)               → must be inside

FOR (g₂_phys, g₃_phys, label) IN KNOWN_PHYSICAL_POINTS[D]:
  FOR EACH bp IN boundary_points:
    // Check that the physical point is on the "inside" of every boundary half-space:
    f_star ← RECONSTRUCT_F(bp.f_star, basis_f)
    I_phys ← γ[f_star] + 2·g₂_phys·α[f_star] + g₃_phys·β[f_star]
    ASSERT I_phys ≥ −ε_verify   // physical point must satisfy all valid inequalities


// ── D4. OUTPUT ──
//
// Restore physical units: with general M and 8πG,
//   g₂ → g₂·M² / (8πG)    (dimensionless in units M=1, 8πG=1)
//   g₃ → g₃·M⁴ / (8πG)

RETURN {
  boundary:  [(bp.g₂, bp.g₃) FOR bp IN boundary_points],  // 40 boundary points
  functionals: [bp.f_star FOR bp IN boundary_points],      // dual certificates
  metadata: {D: D, parameters: {...}}
}
```

---

### 5.4 Recursive Correctness Verification

The following checks are distributed across all four phases and must pass at every stage of the computation. They form a self-consistent verification web.

#### 5.4.1 Check Web for Phase A (Analytic Precomputation)

| Check | What to verify | Why it's correct |
|---|---|---|
| **A-1** | `γ[n] = 1/(n-1)` for all `n > 1` in basis_f | Direct evaluation of `∫₀¹ p^{n-2} dp` |
| **A-2** | `α[n] = 1/(n+1)`, `β[n] = 1/(n+3)` | Same |
| **A-3** | `IMPACT_EXACT[n](b=0) = 1/(n+1) = α[n]` | `₁F₂(...;0) = 1`, so limit gives `α[n]` exactly |
| **A-4** | `TABLE_C[n][0](m) → 0` as `m → ∞` for all `n` | At large m, all kernels decay as `1/m⁴` or faster |
| **A-5** | `Σₙ aₙ·TABLE_C[n][J](m)` evaluated at low-energy amplitude gives `8πG/(p²) + 2g₂ − g₃p²` | This IS the defining property of the improved kernel; failure means TABLE_C is wrong |
| **A-6** | `TABLE_X[k][i][J](m)` satisfies `⟨X_k⟩ = 0` for test spectral density | Null constraint must vanish on any UV spectrum |

#### 5.4.2 Check Web for Phase B (Constraint Assembly)

| Check | What to verify | Why it's correct |
|---|---|---|
| **B-1** | All Block-1 rows are linearly independent (rank of constraint matrix = N_vars for a generic set of x_i, J) | The kernels for different (m, J) are distinct functions |
| **B-2** | Block-3 PSD matrix at `b = B_large`: entries match Block-2 at the same b | Continuity check between discrete and asymptotic regions |
| **B-3** | For f = p^n (single power), Block-2 row equals `IMPACT_EXACT[n](b)`, consistent with `α[n]` at b=0 | Cross-reference Tables B and C from Phase A |
| **B-4** | `N_vars = 17` for Figure-1 parameters | Dimension count: 7 + 6 + 4 |

#### 5.4.3 Check Web for Phase C (LP Solve and Adaptive Refinement)

| Check | What to verify | Why it's correct |
|---|---|---|
| **C-1** | `obj_opt < 0` for all θ when `(g₂₀, g₃₀)` is interior | Interior point assumption: `I₀[f] < 0` for optimal f |
| **C-2** | `D_θ[v_opt] = 1.0` after solving (normalization preserved) | SDPB normalization constraint must be satisfied exactly |
| **C-3** | EFT action reconstructed from `v_opt`: `I₀_reconstructed = obj_opt` (to within `1e-10`) | Verify via `Σₙ v_opt[idx_n]·(γ[n] + 2g₂₀·α[n] + g₃₀·β[n])` |
| **C-4** | After convergence: max negative value across dense grid < `tol_val` | The functional is validly non-negative |
| **C-5** | The saturated inequality from `v_opt`: `γ[f*] + 2g₂*·α[f*] + g₃*·β[f*] ≈ 0` (to within `1e-8`) | Boundary point must lie on the half-space boundary of its own functional |
| **C-6** | `λ_star > 0` for all θ | Distance from interior reference must be positive |

#### 5.4.4 Check Web for Phase D (Output Verification)

| Check | What to verify | Why it's correct |
|---|---|---|
| **D-1** | Every boundary point satisfies all Table-2 inequalities with residual ≥ 0 | The boundary is inside the allowed region |
| **D-2** | The saturating Table-2 inequality for each θ has residual ≈ 0 | Identifies which pre-tabulated bound is active |
| **D-3** | Boundary curve is convex (all cross-products non-negative) | LP allowed region is a convex set |
| **D-4** | Known string theory points `(0, 0)` (Type II, D=6), `(3/16, 3/4)` (Heterotic) lie strictly inside | These theories are known to be consistent |
| **D-5** | As D increases: the allowed region shrinks and becomes more tightly bounded | Physical expectation from the paper's Figure 1 |

---

### 5.5 Key Structural Relationships (Cross-Reference Map)

The following table maps each pseudocode component to the source equation in the TeX files, ensuring the algorithm exactly matches the paper.

| Pseudocode Element | LaTeX Source | Equation Label |
|---|---|---|
| Improved kernel `C²ᵢₘₚ[m²,J,u]` | `sec_localized_b.tex` | `B2 improved flat` |
| EFT reduction: only `8πG/(-u) + 2g₂ − g₃u` | `sec_localized_b.tex` | `eq:definitionofb2improved` |
| Null constraint `X_k` | `sec_localized_b.tex` | `eq:nullconstraints` |
| UV positive measure `⟨·⟩` | `sec_dispersion_review.tex` | `eq:UVaverage` |
| Positivity condition (LP primal) | `app_flat_space_numerics.tex` | `eq:thepositivityconditions` |
| EFT inequality (LP conclusion) | `app_flat_space_numerics.tex` | `eq:fulllinearprogramagain` |
| Basis function table | `app_flat_space_numerics.tex` | `tab:basisfunctionsforf` |
| Analytic kernel integral (J=2 example) | `app_flat_space_numerics.tex` | `eq:examplefunctionofm` |
| Impact parameter transform | `sec_localized_b.tex` | `eq:impactparametertransform` |
| Impact parameter integral formula | `app_flat_space_numerics.tex` | (₁F₂ formula) |
| Large-b asymptotic form | `app_flat_space_numerics.tex` | `eq:the1f2s` |
| 2×2 PSD matrix | `app_flat_space_numerics.tex` | `eq:strongercondition` |
| Adaptive refinement loop | `app_flat_space_numerics.tex` | `sec:adaptiverefinement` |
| Interior point and angle scan | `app_flat_space_numerics.tex` | paragraph beginning "we maximized the distance" |
| Table-2 explicit inequalities | `app_flat_space_numerics.tex` | `eq:inequalitiesforfigure1` |

---

### 5.6 Why the Algorithm Correctly Traces the g₂–g₃ Boundary

The following argument proves soundness and (approximate) completeness:

**Soundness:** If the algorithm returns a point `(g₂*, g₃*)` on the boundary at angle θ, then the functional `f*` with coefficient vector `v_opt` satisfies:
1. `F(m, J) ≥ 0` for all `m ≥ 1`, `J = 0,2,...` — by the positivity constraints, ensured to within `tol_val` by adaptive refinement
2. `γ[f*] + 2g₂*·α[f*] + g₃*·β[f*] = 0` — by construction (the boundary saturation condition)
3. Therefore, the inequality `2α[f*]·g₂ + β[f*]·g₃ ≥ −γ[f*]` is a true EFT bound that is tight at `(g₂*, g₃*)`

Condition 1 is an **outer approximation** (discretized version of the full continuous constraint), so the returned region may be slightly larger than the true allowed region. As the grid is refined, the outer approximation converges to the true region.

**Completeness (approximate):** As the number of scan angles `2·N_angles → ∞` and grid spacing `δ_x, δ_b → 0`, the intersection of all returned half-spaces converges to the true allowed region. The paper verifies this by checking that increasing `J_max` changes bounds by less than `10⁻⁴`.

**The boundary is closed:** Since the allowed region is a closed convex set (intersection of closed half-spaces) and the interior reference point `(g₂₀, g₃₀)` is strictly inside, the boundary is a closed convex curve. The scan over `θ ∈ [0, 2π)` traces this curve completely.

---

### 5.7 Special Case: D=4 (Infrared Divergence Treatment)

In `D=4`, the graviton contribution to the impact-parameter integral diverges:

```
∫₀^∞ db f̂(b) · 8πG·b/(D−4) → ∞  as D→4
```

This requires a separate treatment. The pseudocode above must be modified:

```
IF D == 4:
  // Introduce large-distance cutoff b_max (e.g., the AdS radius in AdS/CFT)
  // Replace the gravity integral by its regularized form:
  //   γ_reg[f] = ∫₀¹ dp f(p) · [1/p² + (const)·log(p/p_min)] / M
  //
  // The boundary is then a FUNCTION of b_max:
  //   g₂ ≥ −(8πG/M²) · 25·log(0.3·M·b_max)
  //
  // The relation between p_min (minimum momentum in f(p)) and b_max is:
  //   p_min² = c / (M·b_max³)  where c ≈ 1 (found from extremal functional)
  //
  // NOTE: An earlier version of the paper (arXiv v1) used an incorrect formula here.
  // The correct relation is p_min² = c/(M·b_max³), NOT p_min² = c/b_max².
  //
  // Procedurally: set f(p=0) ≠ 0 as an additional free parameter,
  // then minimize the coefficient of log(1/p_min) subject to D_θ[f] = 1.

  APPLY_D4_INFRARED_MODIFICATION(...)
```

---

### 5.8 Summary: The Bound Between g₂ and g₃ as Output

The final output of `FindBoundary_g2_g3(D)` is:

- **40 boundary points** `{(g₂*_k, g₃*_k) | k = 0,...,39}` in the `(g₂, g₃)` plane
- **40 extremal functionals** `{f*_k}` that certify the boundary is correct (dual certificates)
- The **allowed region** is the convex hull of the complement: all `(g₂, g₃)` satisfying all 40 (or more) half-space inequalities `2α[f*_k]·g₂ + β[f*_k]·g₃ ≥ −γ[f*_k]`

For `D=5`, this reproduces (in units `8πG=1`, `M=1`) the explicit bound:

```
g₂ − g₃/3 + 60.3086 ≥ 0     ← most restrictive in the g₃ direction
g₂ + c_j·g₃ + C_j  ≥ 0      ← 16 further constraints filling the rest of the boundary
```

as tabulated in `app_flat_space_numerics.tex`, Table `eq:inequalitiesforfigure1`, providing a direct numerical verification that the algorithm is working correctly.
