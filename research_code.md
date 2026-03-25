# Code Documentation Report: `test5.m`

**Subject:** Gravitational EFT Bootstrap — Bounding `g₂` and `g₃` via S-matrix positivity  
**Paper:** Caron-Huot, Mazac, Rastelli, Simmons-Duffin — *Bounding the Space of Gravitational EFT Couplings with S-Matrix Positivity*  
**Solver:** SDPB 3.0.0, via the `SDPB.m` Mathematica package (`WritePmpJson`)  
**Units throughout:** `M = 1`, `8πG = 1` (set in the paper's appendix convention, restored by dimensional analysis)

---

## Table of Contents

1. [Global Parameters: `D`, `Nf`](#1-global-parameters)
2. [Section 1 — Gegenbauer Polynomial: `PJ`, `PJprime1`](#2-section-1--gegenbauer-polynomial)
3. [Section 2 — Improved Kernel: `T1`, `T2`, `T3`, `KernelC`](#3-section-2--improved-kernel)
4. [Section 3 — EFT Action Integrals: `gamma`, `alpha`, `beta`](#4-section-3--eft-action-integrals)
5. [Section 4 — D=5 Difference Basis](#5-section-4--d5-difference-basis)
6. [Section 5 — Null Constraint Kernel: `NhByK`, `XKernel`, `KernelX`](#6-section-5--null-constraint-kernel-and-nhbyk)
7. [Section 6 — Impact Parameter Integrals](#7-section-6--impact-parameter-integrals)
8. [Section 7 — Discretization Parameters and Variable Layout](#8-section-7--discretization-parameters)
9. [Section 8 — Memoized Numerical Evaluation](#9-section-8--memoized-numerical-evaluation)
10. [Section 9 — Constraint Row Assembly](#10-section-9--constraint-row-assembly)
11. [Section 10 — Objective and Normalization Vectors](#11-section-10--objective-and-normalization-vectors)
12. [Section 11 — SDPB Interface: `BuildAndWriteSDP`](#12-section-11--sdpb-interface)
13. [Section 12 — Solution Reconstruction: `ReadSDPBSolution`](#13-section-12--solution-reconstruction)
14. [Section 13 — Adaptive Refinement Loop](#14-section-13--adaptive-refinement-loop)
15. [Section 14 — Boundary Extraction: `ExtractBoundaryPoint`](#15-section-14--boundary-extraction)
16. [Section 15 — Full Boundary Scan: `FindBoundary`](#16-section-15--full-boundary-scan)
17. [Section 16 — Verification Against Paper Table 2](#17-section-16--verification)
18. [End-to-End Data Flow](#18-end-to-end-data-flow)
19. [SDPB Manual Connection Summary](#19-sdpb-manual-connection-summary)

---

## 1. Global Parameters

```mathematica
D  = 5;
Nf = 7;
```

`D` is the spacetime dimension. Every analytic formula in the file — Gegenbauer polynomials, the improved kernel, impact parameter integrals, asymptotic coefficients — depends on `D` algebraically. Setting it as a global before any function is defined means Mathematica substitutes it at call time via its `DownValues` mechanism. Changing `D` to 6, 7, …, 12 re-runs the same algorithm for different dimensions, reproducing the other panels of Figure 1 from the paper.

`Nf = 7` sets the number of basis functions used to expand `f(p)`. For D=5 with the difference basis `pᵏ − p²`, the powers run from `k=3` to `k=Nf+2=9`, giving 7 elements: `{p³−p², p⁴−p², …, p⁹−p²}`. This is the "7.5" entry in the paper's Table 1 notation (the row labelled `{pⁿ C₂^{improved} | n ≤ 7.5}` in `app_flat_space_numerics.tex`).

Combined with `N_h[4] = 6` and `N_h[6] = 4` null-constraint powers (discussed in Section 6 below), this gives **7 + 6 + 4 = 17 total decision variables** — exactly the "17-dimensional space of functionals" cited in the paper's Table 1 for Figure 1.

---

## 2. Section 1 — Gegenbauer Polynomial

```mathematica
PJ[J_, x_] :=
  GegenbauerP[J, (D-3)/2, x] / GegenbauerP[J, (D-3)/2, 1];

PJprime1[J_] := J*(J + D - 3) / (D - 2);
```

### Physics

The unitarity decomposition of the 2→2 scattering amplitude requires partial waves expanded in angular functions on the sphere `S^{D-2}`. These are the Gegenbauer polynomials. The paper's convention (used throughout `sec_localized_b.tex`, including eq. B2 improved flat) is:

```
P_J(x) = ₂F₁(−J, J+D−3; (D−2)/2; (1−x)/2)
```

This specific normalization satisfies `P_J(1) = 1` for all `J`, which is important because `P_J(1) = 1` appears explicitly in the second term of the improved kernel (the T2 term).

### Mathematica implementation

Mathematica's built-in `GegenbauerP[J, λ, x]` computes the standard `C_J^λ(x)` with `λ = (D−3)/2`. This is proportional to the paper's `P_J(x)` but normalized differently: `C_J^λ(1) = Γ(2λ+J)/(Γ(2λ)Γ(J+1)) ≠ 1` in general. Dividing by `GegenbauerP[J,(D-3)/2,1]` enforces the paper's normalization.

**Why this matters:** If the raw `GegenbauerP` were used without dividing by its value at 1, then `T2` (the `P_J(1)` term) would not reduce to a simple `J`-independent expression. The entire three-way split into T1/T2/T3 depends on `P_J(1) = 1`.

### Derivative at `x = 1`

The paper's eq. B2 improved flat contains the term `4u·P_J'(1)/(m⁴−u²)`. The derivative of the normalized Gegenbauer at `x=1` satisfies the closed-form identity:

```
P_J'(1) = J(J+D−3)/(D−2)
```

This is related to the quadratic Casimir of `SO(D−1)`. The key point is that only this single value — `P_J'` evaluated at `x=1`, not the full derivative function — ever appears. So `PJprime1[J]` is stored as a simple algebraic expression in `J` and `D`, requiring no symbolic differentiation inside any integral.

---

## 3. Section 2 — Improved Kernel

```mathematica
T1[n_, J_, m_] := Integrate[
  p^n * (2*m^2 - p^2) * PJ[J, 1 - 2*p^2/m^2] / (m^2 * (m^2 - p^2)^2),
  {p, 0, 1}, Assumptions -> {m > 1, n > 1}];

T2[n_, m_] := Integrate[
  -p^n * p^4 * (4*m^2 - 3*p^2) / (m^6 * (m^2 - p^2)^2),
  {p, 0, 1}, Assumptions -> {m > 1, n > 1}];

T3[n_, J_, m_] := PJprime1[J] * Integrate[
  4 * p^n * p^6 / (m^6 * (m^4 - p^4)),
  {p, 0, 1}, Assumptions -> {m > 1, n > 1}];

KernelC[n_, J_, m_] := T1[n, J, m] + T2[n, m] + T3[n, J, m];
```

### Source: eq. B2 improved flat (`sec_localized_b.tex`)

The improved sum rule kernel is:

```
C²ᵢₘₚ[m², J, u] = (2m²+u)·P_J(1+2u/m²)          u²   (4m²+3u)·P_J(1)     4u·P_J'(1)
                   ───────────────────────────  −  ─── · ─────────────────  + ────────────
                        m²·(m²+u)²                 m⁶       (m²+u)²             m⁴−u²
```

`KernelC[n, J, m]` computes the integral `∫₀¹ dp pⁿ · C²ᵢₘₚ[m², J, −p²]`, where the substitution `u = −p²` is applied throughout (momentum transfer is spacelike).

The formula is split into three named pieces because each has a different mathematical structure that Mathematica's `Integrate` handles separately.

### T1: the J-dependent main term

Under `u = −p²`:
- `(2m²+u) → (2m²−p²)`
- `(m²+u)² → (m²−p²)²`
- `P_J(1+2u/m²) → P_J(1−2p²/m²)` — a Gegenbauer polynomial evaluated at `1−2p²/m²`

The integrand of T1 is `pⁿ·(2m²−p²)·P_J(1−2p²/m²) / [m²·(m²−p²)²]`. Mathematica evaluates this symbolically as a `₂F₁` hypergeometric function of `1/m²`, which matches the example formula in eq. `examplefunctionofm` of `app_flat_space_numerics.tex`.

### T2: the J-independent term

The second group of the kernel contains `P_J(1) = 1` (by our normalization), making it completely independent of `J`. Under `u = −p²`:
- `u²/m⁶ → p⁴/m⁶`
- `(4m²+3u) → (4m²−3p²)`
- `(m²+u)² → (m²−p²)²`

The integrand is `−pⁿ·p⁴·(4m²−3p²) / [m⁶·(m²−p²)²]`. The negative sign comes directly from the minus sign in front of the second group in the kernel formula. T2 produces rational/logarithmic expressions after integration.

### T3: the Casimir-weighted term

The remaining piece of the second group is `4u·P_J'(1)/(m⁴−u²)`. Under `u = −p²`:
- `u² → p⁴`, so `u²/m⁶ → p⁴/m⁶`
- `4u → −4p²`
- `m⁴−u² → m⁴−p⁴`

Assembling the sign: the formula has `−(u²/m⁶)·4u·P_J'(1)/(m⁴−u²)`. Substituting:

```
−(p⁴/m⁶)·(−4p²)·P_J'(1)/(m⁴−p⁴) = +4p⁶·P_J'(1)/[m⁶·(m⁴−p⁴)]
```

**The sign is positive.** This is easy to get wrong by not tracking the double sign flip from `u² → p⁴` (positive) and `4u → −4p²` (negative), whose product is negative, which then cancels the outer minus sign.

`PJprime1[J]` is factored out entirely because it is a scalar at fixed `J`, making the remaining integral `J`-independent. This saves computation when tabulating over many `J` values.

### Assumptions

`Assumptions -> {m > 1, n > 1}` are necessary. `m > 1` tells Mathematica the mass is above threshold, so `(m²−p²) > 0` for `p ∈ [0,1]` — the denominator is non-zero and the integral is well-defined. `n > 1` ensures `∫₀¹ pⁿ/p² dp = ∫₀¹ p^{n-2} dp = 1/(n-1)` converges (the gravity coupling `1/p²` would diverge otherwise).

---

## 4. Section 3 — EFT Action Integrals

```mathematica
gamma[n_] := 1/(n - 1);
alpha[n_] := 1/(n + 1);
beta[n_]  := 1/(n + 3);
```

### Source: eq. `fulllinearprogramagain` (`app_flat_space_numerics.tex`)

The EFT side of the master inequality is (in units `M=1`, `8πG=1`):

```
∫₀¹ dp f(p) · [1/p² + 2g₂ + g₃p²] = γ[f] + 2g₂·α[f] + g₃·β[f]
```

where:
```
γ[f] = ∫₀¹ dp f(p)/p²       (coefficient of the gravity/Newton term)
α[f] = ∫₀¹ dp f(p)           (coefficient of 2g₂)
β[f] = ∫₀¹ dp f(p)·p²       (coefficient of g₃)
```

For a pure power `f(p) = pⁿ`, these are trivial integrals with no computation needed:
```
γ[n] = ∫₀¹ p^{n-2} dp = 1/(n−1)   (valid only for n > 1)
α[n] = ∫₀¹ pⁿ dp = 1/(n+1)
β[n] = ∫₀¹ p^{n+2} dp = 1/(n+3)
```

These three functions are the building blocks of both the objective vector (Section 10) and the boundary extraction formula (Section 14). They are exact rational numbers — no floating point is introduced here.

**Self-check:** `ImpactIntegral[n, 0] = HypergeometricPFQ[{(n+1)/2},{...},0]/(n+1) = 1/(n+1) = alpha[n]`, confirming that the impact parameter integral at zero impact parameter recovers the forward-limit `g₂` action, as it must physically.

---

## 5. Section 4 — D=5 Difference Basis

```mathematica
basisPowers = Range[3, Nf + 2];   (* {3,4,5,6,7,8,9} for Nf=7 *)

gamma_d[k_] := gamma[k] - gamma[2];   (* 1/(k-1) - 1  *)
alpha_d[k_] := alpha[k] - alpha[2];   (* 1/(k+1) - 1/3 *)
beta_d[k_]  := beta[k]  - beta[2];    (* 1/(k+3) - 1/5 *)

KernelC_d[k_, J_, m_] := KernelC[k, J, m] - KernelC[2, J, m];
```

### Why differences are needed in D=5

From `app_flat_space_numerics.tex` (lines 97–98 and Table 1), for each dimension D there is a specific basis for `f(p)`. For D=5 specifically, the basis is `{pᵏ−p² | k=3,4,5,…}` rather than pure powers `pᵏ`.

The reason is explained by the large-`b` asymptotic of the impact parameter integral (eq. `the1f2s`):

```
ImpactIntegral[n,b] ~ A_n/b^{n+1} + B_n·cos(b − π(D−1)/4)/b^{(D-1)/2} + …
```

where `B_n = 2^{(D-3)/2}·Γ((D-2)/2)/√π` is **the same for all values of `n`**. It depends on `D` but not on the exponent `n`.

For a sum `Σ aₙ pⁿ`, the oscillatory part of the impact parameter transform is `(Σ aₙ)·B·cos(b−…)/b^{(D-1)/2}`. For D=5, `(D−1)/2 = 2`, so this oscillatory term decays as `1/b²`. The non-oscillatory term for a single power `pⁿ` decays as `1/b^{n+1}`. But for D=5 and any `n ≥ 2`, we have `n+1 ≥ 3 > 2`, meaning the **oscillatory term dominates at large b**. A function that has a leading `cos(b)/b²` term cannot be globally non-negative for all `b`, because cosine alternates in sign.

The paper's solution: use differences `pᵏ − p²`. The oscillatory `B`-coefficient satisfies:
```
B_{pᵏ-p²} = B_k − B_2 = 0
```
because `B_n` is n-independent. The leading oscillatory term exactly cancels, leaving only the non-oscillatory term `(A_k − A_2)/b^{k+1}`, which can be made positive.

For odd `D ≥ 7`, pure powers `p², p³, p⁴, …` work because `(D−3)/2 ≥ 2`, meaning there exist powers `n ≤ (D-3)/2` whose non-oscillatory term dominates. For even `D ≥ 6`, half-integer powers `p^{3/2}, p^{5/2}, …` are used for similar reasons. D=5 is the exceptional case requiring differences.

### Indexing

`basisPowers = {3,4,5,6,7,8,9}` for `Nf=7`. The lower bound is 3 because `k=2` would give `p² − p² = 0`. The upper bound is `Nf+2 = 9`.

The difference action integrals are differences of the pure-power versions:
```
gamma_d[k] = 1/(k-1) − 1       (e.g. for k=3: 1/2 − 1 = −1/2)
alpha_d[k] = 1/(k+1) − 1/3     (e.g. for k=3: 1/4 − 1/3 = −1/12)
beta_d[k]  = 1/(k+3) − 1/5     (e.g. for k=3: 1/6 − 1/5 = −1/30)
```

These small numbers control the shape of the allowed region in the `(g₂, g₃)` plane. The objective and normalization vectors are assembled from them.

---

## 6. Section 5 — Null Constraint Kernel and `NhByK`

```mathematica
NhByK = <|4 -> 6, 6 -> 4|>;
```

### What `NhByK` is

`NhByK` is a **Mathematica Association** (a key-value dictionary). The syntax `<| key -> value, … |>` creates it. This specific Association encodes the number of basis powers used to expand each null constraint function `h_k(p)`:

- `NhByK[4] = 6`: the function `h_4(p)` is expanded in `{p⁰, p¹, p², p³, p⁴, p⁵}` — **six powers** (indices `i = 0, 1, 2, 3, 4, 5`)
- `NhByK[6] = 4`: the function `h_6(p)` is expanded in `{p⁰, p¹, p², p³}` — **four powers** (indices `i = 0, 1, 2, 3`)

### Where it comes from in the paper

The paper's Table 1 (`app_flat_space_numerics.tex`) lists the full set of functionals used for Figure 1:

```
{pⁿ C₂^{improved} | n ≤ 7.5}  ∪  {pⁱ X₄ | i = 0,1,2,3,4,5}  ∪  {pⁱ X₆ | i = 0,1,2,3}
```

The second group gives 6 null-constraint variables for `k=4`, and the third group gives 4 for `k=6`. So `NhByK` is a direct transcription of those two rows of Table 1: `4 → 6` and `6 → 4`.

### Why the null constraints improve the bounds

The null constraints `⟨X_{k,u}[m², J]⟩ = 0` (eq. nullconstraints in `sec_localized_b.tex`) are exact consequences of crossing symmetry. They hold for any physical UV completion. Because they equal zero on the UV side, multiplying by any function `h_k(p)` and integrating adds zero to the UV side of the master inequality — but it can change the EFT side by contributing terms involving `h_k`. Including them allows the LP to use more degrees of freedom to sharpen the allowed region without violating any physical principle.

Without null constraints (i.e., if `NhByK` were empty), the optimization would only use 7 variables and produce weaker bounds. The 10 extra variables from the null constraints (6 for `k=4` and 4 for `k=6`) tighten the boundary of the allowed region in Figure 1.

### How `NhByK` is used throughout the code

The Association is read in four places:

**1. Total variable count (Section 7):**
```mathematica
Nvars = Nf + Total[Values[NhByK]]   (* = 7 + 6 + 4 = 17 *)
```
`Total[Values[NhByK]]` sums the values `{6, 4}` to give 10, the total number of null-constraint variables.

**2. Index helpers (Section 7):**
```mathematica
h4Idx[i_] := Nf + 1 + i    (* i=0..5 → positions 8..13 in v *)
h6Idx[i_] := Nf + 7 + i    (* i=0..3 → positions 14..17 in v *)
```
These map `h_k` power indices to positions in the decision variable vector `v`.

**3. Block 1 row assembly (Section 9):**
```mathematica
Do[row[[h4Idx[i]]] = KernelXNum[4, i, J, mval], {i, 0, NhByK[4]-1}];
Do[row[[h6Idx[i]]] = KernelXNum[6, i, J, mval], {i, 0, NhByK[6]-1}];
```
The loops run `i = 0` to `NhByK[4]-1 = 5` and `i = 0` to `NhByK[6]-1 = 3`, filling the constraint row's null-constraint entries.

**4. Functional evaluation (Section 13):**
```mathematica
Sum[v[[h4Idx[i]]] * KernelXNum[4, i, J, m], {i, 0, NhByK[4]-1}]
Sum[v[[h6Idx[i]]] * KernelXNum[6, i, J, m], {i, 0, NhByK[6]-1}]
```
When validating the solution, the null-constraint contributions are re-evaluated using the same loop bounds from `NhByK`.

**5. Block 3 PSD polynomial (Section 11):**
```mathematica
ConstantArray[{0}, NhByK[4] + NhByK[6]]
```
The h-variables do not appear in the large-`b` asymptotic (they are subleading as `m→∞`), so their polynomial contribution is set to zero. The total count `NhByK[4] + NhByK[6] = 10` determines how many zero entries to append.

### The null constraint formula

The full formula (eq. nullconstraints in `sec_localized_b.tex`) is:

```
X_k[m², J, u] = [(2m²+u)·m²·P_J(1+2u/m²)] / (u·m²·(m²+u))^{1+k/2}

               − Res_{u'=0} [ (2m²+u')(m²-u')(m²+2u') · m²·P_J(1+2u'/m²) ]
                             / [ m²·(u-u')·u'·(m²-u)·(m²+u')·(m²+u+u') · (u'·m²·(m²+u'))^{k/2} ]
```

`XKernelFirstTerm` implements the first line. `XKernelResidueIntegrand` implements the integrand of the residue, and `Residue[..., {up, 0}]` extracts the coefficient of `1/up` in its Laurent expansion at `up=0`. This residue subtracts all EFT contact contributions, leaving a purely UV quantity that equals zero on any physical spectrum.

The pole order at `u'=0` is `1 + k/2` (three for `k=4`, four for `k=6`), matching the paper's description.

`KernelX[k, i, J, m]` then integrates `pⁱ · X_k[m², J, −p²]` over `p ∈ [0,1]`. This gives the Table B entry in the pseudocode: the numerical coefficient in the constraint matrix for null-constraint variable `b_{k,i}` at mass point `m` and spin `J`.

---

## 7. Section 6 — Impact Parameter Integrals

```mathematica
ImpactIntegral[n_, b_] :=
  HypergeometricPFQ[{(n+1)/2}, {(D-2)/2, (n+3)/2}, -b^2/4] / (n+1);

ImpactIntegral_d[k_, b_] :=
  ImpactIntegral[k, b] - ImpactIntegral[2, b];

ImpactCoeffA[n_]   := 2^n * Gamma[(D-2)/2] * Gamma[(n+1)/2] / Gamma[(D-n-3)/2];
ImpactCoeffB[n_]   := 2^((D-3)/2) * Gamma[(D-2)/2] / Sqrt[Pi];
ImpactCoeffA_d[k_] := ImpactCoeffA[k] - ImpactCoeffA[2];
```

### Derivation of the exact formula

The impact parameter positivity condition comes from taking the `m → ∞` limit of the kernel constraint with `b = 2J/m` held fixed. The Gegenbauer-to-Bessel limit (eq. GegenbauerToBessel in `sec_localized_b.tex`) says:

```
lim_{m→∞} P_{mb/2}(1 − 2p²/m²) = Γ((D-2)/2)·J_{(D-4)/2}(bp) / (bp/2)^{(D-4)/2}
```

The positivity condition on `f(p)` therefore requires:

```
Γ((D-2)/2) · ∫₀¹ dp f(p)·J_{(D-4)/2}(bp)/(bp/2)^{(D-4)/2} ≥ 0   for all b ≥ 0
```

For `f(p) = pⁿ`, this Bessel integral has the closed form (eq. the1f2s in `app_flat_space_numerics.tex`):

```
₁F₂((n+1)/2 ; (D-2)/2, (n+3)/2 ; −b²/4) / (n+1)
```

This is exactly `ImpactIntegral[n, b]`. Mathematica's `HypergeometricPFQ` evaluates this to arbitrary precision.

The null constraints `X_k` do not appear here because they are subleading in the `m → ∞` limit, so `h_k(p)` drops out — this is why Block 2 rows contain only f-variable entries.

### Asymptotic coefficients

`ImpactCoeffA[n]` is the coefficient of the non-oscillatory `1/b^{n+1}` term in the large-`b` expansion, and `ImpactCoeffB[n]` is the coefficient of the oscillatory `cos(b−π(D-1)/4)/b^{(D-1)/2}` term. As discussed in Section 5, `B_n` is `n`-independent, so `ImpactCoeffA_d[k] = ImpactCoeffA[k] - ImpactCoeffA[2]` is the non-oscillatory coefficient of the difference basis function `pᵏ − p²` after the B-cancellation. This appears in the Block 3 polynomial matrix construction.

---

## 8. Section 7 — Discretization Parameters

```mathematica
Jmax     = 42;
Jlist    = Range[0, Jmax, 2];         (* {0,2,4,...,42}, 22 values *)
deltax   = 1/400;
xgrid0   = N @ Range[0, 1-deltax, deltax];
epsilonb = 1/250;
deltab   = 1/32;
Blarge   = 40;
bgrid0   = N @ Range[epsilonb, Blarge-deltab, deltab];

Nvars  = Nf + Total[Values[NhByK]];   (* = 17 *)

fIdx[ki_]  := Position[basisPowers, ki][[1, 1]];
h4Idx[i_]  := Nf + 1 + i;
h6Idx[i_]  := Nf + 7 + i;
```

### Source: Table 1 of `app_flat_space_numerics.tex`

Every parameter here is taken from the paper's Table 1, column "Figure 1":

| Code symbol | Paper symbol | Value | Meaning |
|-------------|--------------|-------|---------|
| `Jmax` | `J_max` | 42 | Highest spin included in finite sum |
| `deltax` | `δ_x` | 1/400 | Spacing of the mass grid in `x = 1 − 1/m²` |
| `epsilonb` | `ε_b` | 1/250 | Smallest impact parameter sampled |
| `deltab` | `δ_b` | 1/32 | Impact parameter grid spacing |
| `Blarge` | `B` | 40 | Cutoff above which the large-`b` asymptotic PSD matrix is used |
| `prec` | `--precision` | 768 | SDPB arithmetic precision in bits |

### The `x` reparameterization

The paper discretizes `x = 1 − 1/m²` rather than `m` directly. As `m` ranges from 1 (threshold) to `∞`, `x` ranges over `[0, 1)`. The spacing `δ_x = 1/400` gives 400 initial grid points: `{0, 1/400, 2/400, …, 399/400}`. The conversion `mFromX[x] = 1/√(1−x)` is used everywhere mass values are needed.

### Only even spins

`Jlist = Range[0, 42, 2] = {0, 2, 4, …, 42}` contains only even integers because the 2→2 scattering of identical scalars has `s ↔ u` crossing symmetry, which forces the partial-wave expansion to contain only even-spin contributions.

### Variable layout in `v`

The 17-dimensional decision variable vector `v` is laid out as:
```
v[1]  = a₃   (coefficient of p³−p²)
v[2]  = a₄   (coefficient of p⁴−p²)
  ⋮
v[7]  = a₉   (coefficient of p⁹−p²)
v[8]  = b_{4,0}   (coefficient of p⁰ in h₄)
v[9]  = b_{4,1}
  ⋮
v[13] = b_{4,5}
v[14] = b_{6,0}   (coefficient of p⁰ in h₆)
  ⋮
v[17] = b_{6,3}
```

`fIdx[ki]` maps the power `ki` in `{3,…,9}` to its position `1,…,7` using `Position`. `h4Idx[i] = 8+i` and `h6Idx[i] = 14+i` give simple arithmetic formulas.

---

## 9. Section 8 — Memoized Numerical Evaluation

```mathematica
mFromX[x_] := 1 / Sqrt[1 - x];

KernelC_dNum[ki_, J_Integer, m_?NumericQ] :=
  KernelC_dNum[ki, J, m] = N[KernelC_d[ki, J, m], 50];

KernelXNum[k_, i_Integer, J_Integer, m_?NumericQ] :=
  KernelXNum[k, i, J, m] = N[KernelX[k, i, J, m], 50];

ImpactNum[ki_, b_?NumericQ] :=
  ImpactNum[ki, b] = N[ImpactIntegral_d[ki, b], 50];
```

### Mathematica memoization

The pattern `f[args_] := f[args] = expr` is the standard Mathematica memoization idiom. The first time `f` is called with specific numeric arguments, it evaluates `expr`, then stores the result as a new `DownValue` `f[specific_args] = result`. All subsequent calls with the same arguments match the stored value directly and skip the evaluation.

The pattern restrictions `J_Integer` and `m_?NumericQ` ensure memoization only triggers for fully numeric inputs — not symbolic ones — preventing premature evaluation during function definition.

### Precision: 50 decimal digits

`N[expr, 50]` evaluates to 50 decimal digits. SDPB runs at 768 bits `≈ 231` decimal digits internally, but the constraint matrix entries fed into SDPB only need to be accurate at the LP level. 50 digits is more than sufficient to avoid numerical noise in the linear program while keeping Mathematica evaluation times manageable.

---

## 10. Section 9 — Constraint Row Assembly

```mathematica
Block1Row[x_, J_] := Module[{row = ConstantArray[0, Nvars], mval = mFromX[x]},
  Do[row[[fIdx[ki]]]  = KernelC_dNum[ki, J, mval], {ki, basisPowers}];
  Do[row[[h4Idx[i]]]  = KernelXNum[4, i, J, mval], {i, 0, NhByK[4]-1}];
  Do[row[[h6Idx[i]]]  = KernelXNum[6, i, J, mval], {i, 0, NhByK[6]-1}];
  row];

Block2Row[b_] := Module[{row = ConstantArray[0, Nvars]},
  Do[row[[fIdx[ki]]]  = ImpactNum[ki, b], {ki, basisPowers}];
  row];

AssembleConstraints[xgrid_, bgrid_] := Module[{rows = {}},
  Do[Do[AppendTo[rows, Block1Row[x, J]], {x, xgrid}], {J, Jlist}];
  Do[AppendTo[rows, Block2Row[b]], {b, bgrid}];
  rows];
```

### The constraint structure

The master LP (eq. `thepositivityconditions` in `app_flat_space_numerics.tex`) requires:

```
F(m, J) ≡ ∫₀¹ dp f(p)·C²ᵢₘₚ[m²,J,−p²]
          + Σ_{k=4,6} ∫₀¹ dp h_k(p)·X_k[m²,J,−p²]  ≥  0
```

for all `m ≥ 1`, `J = 0,2,4,…`. Expanding `f(p) = Σ aₖ(pᵏ−p²)` and `h_k(p) = Σᵢ b_{k,i}·pⁱ`, this becomes a linear inequality in the 17-dimensional vector `v`:

```
row · v ≥ 0
```

where `row[n]` holds the integrated kernel value for variable `n` at a specific `(m, J)` point.

**Block 1 rows** enforce this for each `(x, J)` pair. The 17 entries are: 7 kernel values for the f-variables, 6 null-constraint kernel values for `h₄`, and 4 for `h₆`.

**Block 2 rows** enforce `Σ aₖ·ImpactIntegral_d[k,b] ≥ 0` for each sampled impact parameter `b ≤ Blarge`. The `h_k` entries are zero because the null constraints are subleading in the `m → ∞` limit (stated explicitly in `app_flat_space_numerics.tex` line 48).

The initial constraint count is `400 × 22 (Block 1) + ≈1279 (Block 2) ≈ 10079` rows.

---

## 11. Section 10 — Objective and Normalization Vectors

```mathematica
g20 = -5;
g30 = -15;

ObjVector[g20_, g30_] := Module[{v = ConstantArray[0, Nvars]},
  Do[v[[fIdx[ki]]] = N[gamma_d[ki] + 2*g20*alpha_d[ki] + g30*beta_d[ki], 50],
     {ki, basisPowers}]; v];

NormVector[theta_] := Module[{v = ConstantArray[0, Nvars]},
  Do[v[[fIdx[ki]]] = N[2*Cos[theta]*alpha_d[ki] + Sin[theta]*beta_d[ki], 50],
     {ki, basisPowers}]; v];
```

### The angle-scan LP formulation

The paper (`app_flat_space_numerics.tex`, first paragraph) describes maximizing the distance from a chosen interior point `(g₂₀, g₃₀)` along rays of constant angle `θ`. This is implemented as:

```
minimize   I₀[f] = γ[f] + 2g₂₀·α[f] + g₃₀·β[f]
subject to D_θ[f] = 2cosθ·α[f] + sinθ·β[f] = 1
           F(m, J) ≥ 0  for all (m, J) on the grid
```

`ObjVector` computes the objective coefficients `obj[k] = gamma_d[k] + 2·g20·alpha_d[k] + g30·beta_d[k]` for each basis power `k`. `NormVector` computes the normalization coefficients `norm[k] = 2·cos(θ)·alpha_d[k] + sin(θ)·beta_d[k]`. Both vectors are zero for the `h_k` variables (positions 8–17), because the null constraints have zero EFT action — `⟨X_k⟩ = 0` identically.

The interior reference point `(g₂₀, g₃₀) = (−5, −15)` is chosen to lie inside the expected allowed region near its tip, so that the scan over all 40 angles traces the full closed boundary curve.

---

## 12. Section 11 — SDPB Interface

```mathematica
<<"../SDPB.m";
prec = 768;

BuildAndWriteSDP[theta_, outFile_String] := Module[...
  sdpObj = SDP[objFull, normFull, Join[scalBlocks, {psdBlock}]];
  WritePmpJson[outFile, sdpObj, prec];
];
```

### SDPB 3.0 PMP convention (Manual Section 3.2)

The SDPB Mathematica package provides `WritePmpJson[file, SDP[obj, norm, posMatrices], prec]`. The SDP object corresponds to the Polynomial Matrix Program (PMP) format (SDPB manual, eq. 3.1):

```
maximize  a·z   subject to:
  Σ_n z_n · W^n_j(x) ⪰ 0  for all x ≥ 0, j = 1…J
  n·z = 1
```

Our mapping:
- `z = (z₀, z₁, …, z₁₇)` where `z₀` is eliminated by the normalization and `z₁…z₁₇ = v[1]…v[17]`
- `objFull = Prepend[-ObjVector[g20,g30], 0]` — negated because SDPB **maximizes** while we want to **minimize** `I₀`
- `normFull = Prepend[NormVector[theta], 0]` — enforces `D_θ[f] = 1`

### Block 1+2: scalar (1×1) PSD blocks

Each constraint row `row·v ≥ 0` becomes one `PositiveMatrixWithPrefactor` with a 1×1 polynomial matrix. Since the constraint has no dependence on the continuous variable `x`, each matrix entry polynomial is degree-0 (a constant). In SDPB's JSON representation, the polynomial vector for the `(1,1)` entry is:

```mathematica
PositiveMatrixWithPrefactor[<|
  "polynomials" -> {{{ Prepend[Table[{rows[[r,n]]}, {n,1,Nvars}], {0}] }}}
|>]
```

The nesting is: `polynomials` → list of polynomial matrices → matrix rows → matrix columns → polynomial vector (one polynomial per decision variable). `{0}` is the degree-0 polynomial `0` for `z₀`; `{rows[[r,n]]}` is the degree-0 polynomial with coefficient `rows[[r,n]]` for variable `z_n`.

### Block 3: 2×2 polynomial PSD matrix

For `b ≥ Blarge`, the paper's eq. `strongercondition` requires:

```
[[A(b)+B(b), C(b)], [C(b), A(b)−B(b)]] ⪰ 0
```

For D=5 with the difference basis, `B(b) = 0` exactly, so this reduces to `A(b)·I₂ ⪰ 0`. After substituting `z = 1/b` and multiplying by `b² = z⁻²` to clear negative powers, the entry `A(b)·b²` becomes:

```
A(b)·b² = Σₖ aₖ · ImpactCoeffA_d[k] · z^{k-1}
```

This is a polynomial in `z` of degree up to `9−1 = 8`. For decision variable `a_k`, the polynomial is a list of 9 coefficients with `ImpactCoeffA_d[k]` at position `k−1` and zeros elsewhere. The `h_k` variables get the zero polynomial `{0}` because they do not appear in the large-`b` limit. Both diagonal entries carry the same polynomial; off-diagonal entries are zero.

This 2×2 polynomial matrix PSD block encodes the constraint that `A(b) ≥ 0` for all `b ≥ Blarge`, replacing infinitely many scalar inequalities with one finite SDP block that SDPB handles natively.

---

## 13. Section 12 — Solution Reconstruction

```mathematica
ReadSDPBSolution[yFile_String, normVec_] := Module[
  {yLines, yVals, k0, zFull, vFull},
  yLines = Import[yFile, "Lines"];
  yVals  = ToExpression /@ Select[yLines, StringMatchQ[#, NumberString] &];
  k0    = First @ Ordering[-Abs[normVec], 1];
  zFull = ConstantArray[0, Length[normVec]];
  Do[Module[{idx = If[i < k0, i, i+1]}, zFull[[idx]] = yVals[[i]]],
     {i, Length[yVals]}];
  zFull[[k0]] = (1 - normVec.zFull) / normVec[[k0]];
  Rest[zFull]
];
```

### SDPB output: `y.txt` (Manual Section 5.3)

After running with `--writeSolution=y,z`, SDPB writes the reduced solution vector `y` to `out_dir/y.txt`. SDPB eliminates the component `z_{k₀}` corresponding to the largest `|n_{k₀}|` in the normalization vector, for numerical stability (SDPB manual, Section 3, paragraph after eq. 3.1). The remaining components form `y`, a vector of length `Nvars` (= 17 here).

`k0 = First @ Ordering[-Abs[normVec], 1]` finds the index of the largest `|n_k|`. The reconstruction loop inserts `y[i]` into position `i` if `i < k0`, or position `i+1` if `i ≥ k0`, effectively inserting a gap at position `k0`. Then `zFull[[k0]] = (1 − n·z_rest)/n_{k0}` recovers the missing component from the normalization equation `n·z = 1`. Finally `Rest[zFull]` drops `z₀` and returns the 17-dimensional `v`.

---

## 14. Section 13 — Adaptive Refinement Loop

```mathematica
AdaptiveRefinement[theta_] := Module[...
  While[!converged && round < 20,
    BuildAndWriteSDP[theta, "pmp_theta.json"];
    Run["pmp2sdp ..."];  Run["sdpb ..."];
    v      = ReadSDPBSolution[...];
    {negX, negB} = ValidateFunctional[v, xgrid, bgrid];
    (* check convergence, refine grids *)
  ];
  {v, objVal}];
```

### Why adaptive refinement is needed

Discretizing `x` with 400 initial points weakens the LP: the solver only enforces `F(m_i, J) ≥ 0` at discrete points, not for all `m`. As the paper explains (`app_flat_space_numerics.tex`, Section adaptiverefinement), the solver typically returns a functional that touches zero at pairs of neighboring grid points and dips slightly negative between them.

The outer-approximation method: locate negative regions, add constraints there, re-solve. Each round the grid grows, the LP tightens, and negative regions shrink by a factor of `Nrefine = 10`.

### `ValidateFunctional`

Evaluates `EvalKernelF[v, J, x]` and `EvalImpact[v, b]` on grids 10× finer than the training grid by constructing `verX` and `verB` from dense sub-intervals around each current grid point. A sign change from non-negative to negative between consecutive `verX` values identifies a negative region `{J, x₁, x₂, value}`.

### `QuadMin`

```mathematica
QuadMin[x1_, f1_, xm_, fm_, x2_, f2_] := Module[{h = xm-x1, A, B},
  A = (f1 - 2*fm + f2) / (2*h^2);
  B = (f2 - f1) / (2*h);
  If[A > 0, Clip[xm - B/(2*A), {x1+1e-10, x2-1e-10}], xm]];
```

Fits a parabola through three equally-spaced points `(x₁,f₁)`, `(xₘ,fₘ)`, `(x₂,f₂)` and returns the minimizer. `A = (f₁ − 2fₘ + f₂)/(2h²)` is the second-difference approximation to `f''`. The minimizer is `x* = xₘ − B/(2A)`. If `A ≤ 0` (concave — unusual but possible), it falls back to the midpoint. This matches the paper's description exactly.

### Refinement grid

For a negative region `(x₁, x₂)` with estimated minimum at `x*`, new points are:
```
{x* − Nrefine·s, …, x*, …, x* + Nrefine·s}   where s = |x₂−x₁|/Nrefine
```
With `Nrefine = 10`, this adds 21 points, reducing the negative interval width by a factor of 10 in the next round.

### Convergence criteria

Convergence is declared when **either**:
- All negative regions have width `< tolNeg = 1e-6` — the functional is essentially everywhere non-negative
- All negative values have `|value| < tolVal = 1e-9` — any remaining negativity is within solver tolerance

The `round < 20` guard prevents infinite loops. The paper notes "typically only a few refinement steps are needed."

---

## 15. Section 14 — Boundary Extraction

```mathematica
ExtractBoundaryPoint[theta_, objOpt_] := Module[{lambda},
  lambda = -objOpt;
  {g20 + lambda * Cos[theta], g30 + lambda * Sin[theta]}];
```

### Derivation

The LP minimizes `I₀[f] = γ[f] + 2g₂₀·α[f] + g₃₀·β[f]` subject to `D_θ[f] = 1`, yielding optimal value `obj_opt`. The optimal functional `f*` defines a supporting half-space boundary:

```
2α[f*]·g₂ + β[f*]·g₃ = −γ[f*]
```

The ray from `(g₂₀, g₃₀)` in direction `(cosθ, sinθ)` intersects this boundary at distance `λ*` where:

```
I₀[f*] + λ*·D_θ[f*] = 0   →   obj_opt + λ*·1 = 0   →   λ* = −obj_opt
```

Since `(g₂₀, g₃₀)` is inside the allowed region, the optimal `f*` has `obj_opt < 0`, so `λ* > 0` — the boundary is strictly outside the reference point in direction `θ`.

**Verification:** substituting back, `2α[f*]·(g₂₀ + λ*cosθ) + β[f*]·(g₃₀ + λ*sinθ) = I₀[f*] − γ[f*] + λ*·1 = obj_opt − γ[f*] + (−obj_opt) = −γ[f*]`. The boundary point exactly saturates the half-space inequality of `f*`. ✓

---

## 16. Section 15 — Full Boundary Scan

```mathematica
Nangles   = 20;
thetaList = Table[k * Pi / Nangles, {k, 0, 2*Nangles-1}];

FindBoundary[] := Module[{pts = {}},
  Do[
    {v, objOpt} = AdaptiveRefinement[theta];
    bpt = ExtractBoundaryPoint[theta, objOpt];
    AppendTo[pts, bpt],
    {theta, thetaList}
  ]; pts];
```

`Nangles = 20` gives `2·Nangles = 40` angles uniformly distributed in `[0, 2π)`, matching the paper's description: "scanned over angles `θ ∈ {0, π/20, …, 39π/20}`." Each angle produces one boundary point; together they trace the full closed boundary curve of the allowed region. The result is a list of 40 `{g₂*, g₃*}` pairs that can be plotted with `ListPlot[pts, Joined→True]`.

---

## 17. Section 16 — Verification

```mathematica
D5Inequalities = {
  {1, -1/3, 60.3086},
  {1, 0.0647867, 9.64034},
  ... (17 rows total) };

CheckD5Point[g2_, g3_] :=
  AllTrue[D5Inequalities, #[[1]]*g2 + #[[2]]*g3 + #[[3]] >= 0 &];
```

### Source

These 17 inequalities are taken verbatim from the paper's Table 2 (eq. `inequalitiesforfigure1` in `app_flat_space_numerics.tex`), D=5 column. Each row `{c₁, c₂, c₃}` encodes `c₁·g₂ + c₂·g₃ + c₃ ≥ 0`. These are the actual numerical output of the paper's computation, listed there as the final result.

After running `FindBoundary[]`, verifying `AllTrue[boundaryPoints, CheckD5Point[#[[1]],#[[2]]]&]` checks that every computed boundary point satisfies all known D=5 inequalities. This is a necessary consistency check: boundary points lie on or inside the allowed region, so they must satisfy every valid inequality.

---

## 18. End-to-End Data Flow

```
GLOBAL SETUP
  D=5, Nf=7
      │
      ├── Sec 1: PJ, PJprime1
      │     Gegenbauer P_J(x) normalized to P_J(1)=1; P'_J(1)=J(J+D-3)/(D-2)
      │
      ├── Sec 2: T1, T2, T3 → KernelC
      │     ∫₀¹ pⁿ·C²ᵢₘₚ[m²,J,−p²] dp  (analytic, via Integrate[])
      │
      ├── Sec 3: gamma, alpha, beta
      │     ∫₀¹ pⁿ·{1/p², 1, p²} dp  (trivial fractions 1/(n-1), 1/(n+1), 1/(n+3))
      │
      ├── Sec 4: basisPowers, gamma_d, alpha_d, beta_d, KernelC_d
      │     D=5 difference basis pᵏ−p²; difference action integrals
      │
      ├── Sec 5: NhByK={4→6, 6→4}, XKernel, KernelX
      │     Null constraint X_k kernel; ∫₀¹ pⁱ·X_k dp  (analytic)
      │
      ├── Sec 6: ImpactIntegral, ImpactCoeffA_d
      │     ₁F₂ impact parameter integrals; asymptotic B-cancellation
      │
      ├── Sec 7: grids, Nvars=17, fIdx/h4Idx/h6Idx
      │     All discretization parameters; variable index layout
      │
      └── Sec 8: KernelC_dNum, KernelXNum, ImpactNum
            Memoized 50-digit numerical evaluation of analytic results

FOR EACH θ ∈ {0, π/20, …, 39π/20}:   [Sec 15: FindBoundary]
      │
      └── AdaptiveRefinement[θ]:        [Sec 13]
            │
            ├── LOOP (up to 20 rounds):
            │     │
            │     ├── Sec 9: AssembleConstraints
            │     │     Build constraint matrix rows (Block1Row + Block2Row)
            │     │     Block1: 400×22=8800 rows, each 17-entry vector
            │     │     Block2: ~1279 rows, f-variables only
            │     │
            │     ├── Sec 10: ObjVector, NormVector
            │     │     17-dim objective (negated) and normalization vectors
            │     │
            │     ├── Sec 11: BuildAndWriteSDP → pmp.json
            │     │     Block1+2 → scalar 1×1 PSD blocks (constant polynomials)
            │     │     Block3  → 2×2 polynomial PSD matrix (large-b asymptotics)
            │     │     WritePmpJson[file, SDP[obj,norm,blocks], 768]
            │     │
            │     ├── SHELL: pmp2sdp → sdp_tmp/
            │     │     Converts JSON to SDPB's internal binary format
            │     │
            │     ├── SHELL: sdpb → sdpb_out/
            │     │     768-bit arithmetic LP solve; writes y.txt, out.txt
            │     │
            │     ├── Sec 12: ReadSDPBSolution
            │     │     Reconstruct 17-dim v from y.txt + normalization n·z=1
            │     │
            │     ├── ValidateFunctional (10× finer grid)
            │     │     Find sign-change intervals negX={J,x₁,x₂,val}, negB={b₁,b₂,val}
            │     │
            │     ├── IF converged (width<1e-6 OR depth<1e-9): BREAK
            │     │
            │     └── QuadMin → add refinement points to xgrid and bgrid
            │
            └── RETURN (v, objOpt)

      └── Sec 14: ExtractBoundaryPoint
            λ* = −objOpt
            (g₂*, g₃*) = (g₂₀ + λ*cosθ, g₃₀ + λ*sinθ)

Sec 16: CheckD5Point
      Verify all 40 boundary points against paper's 17 known inequalities
```

---

## 19. SDPB Manual Connection Summary

| Code element | SDPB concept | Manual location |
|---|---|---|
| `WritePmpJson[file, SDP[...], prec]` | PMP JSON output | Sec. 3.2, Listing 2 |
| `SDP[objFull, normFull, posMatrices]` | PMP problem object | Sec. 3.2 |
| `objFull` | Objective vector `a` | Sec. 3.2, `objective` field |
| `normFull` | Normalization vector `n` (with `n·z=1`) | Sec. 3.2, `normalization` field |
| `PositiveMatrixWithPrefactor[<\|"polynomials"->...\|>]` | Polynomial matrix constraint | Sec. 3.1, eq. 2.1 |
| `scalBlocks` (degree-0, 1×1) | Scalar linear constraints | Sec. 3.1, constant case |
| `psdBlock` (degree-8, 2×2) | Polynomial matrix PSD block | Sec. 3.1, matrix case |
| `Prepend[-ObjVector[...], 0]` | Negated objective (SDPB maximizes) | Sec. 3.1 |
| `prec = 768` | `--precision=768` | Sec. 5.1, Table 1 of paper |
| `Run["pmp2sdp ..."]` | Conversion to SDPB format | Sec. 3.2 |
| `Run["sdpb -s sdp_tmp -o out --writeSolution=y,z"]` | LP solve + solution output | Sec. 5.1, 5.3 |
| `ReadSDPBSolution` reads `y.txt` | Dual solution vector `y` | Sec. 5.3 |
| Reconstruction from `n·z=1` | Elimination of `z_{k₀}` | Sec. 3 (after eq. 3.1) |
| `StringContainsQ[#,"primalObjective"]` | Primal objective from `out.txt` | Sec. 5.3, Listing 7 |
