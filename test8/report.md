# Technical Report: `generate.m` and the Maximal SUSY PMP Problem
## Correspondence between Appendix B of arXiv:2510.07991v1 and the Code

---

## 1. Overview

The file `generate.m` is a Wolfram Mathematica script that constructs a **Polynomial Matrix Program (PMP)** — the input format consumed by the semidefinite program solver SDPB — for the gravitational (maximal SUSY, −2 subtractions) S-matrix bootstrap problem described in **Appendix B** of "The Rise of Linear Trajectories" (Huang, Ricossa, Riva, Tsai, arXiv:2510.07991v1). The script writes its output to `out.json`, which SDPB then optimizes.

The PMP encodes a search for the functional vector **α** that maximizes the coupling of a chosen higher-spin resonance, subject to the positivity constraints imposed by unitarity and crossing symmetry, as formulated in Eqs. (B3)–(B5) of the paper. Below, each component of the code is traced carefully back to the corresponding equation or concept in the appendix.

---

## 2. Physical and Mathematical Setup (Paper Side)

### 2.1 The −2 Subtraction Dispersion Relation and the Graviton Pole

For maximal SUSY (32 supercharges, N=2 in 10D), the amplitude satisfies **−2 subtracted dispersion relations**. The graviton pole, which is unavoidable in this case, appears in the dispersion relation as (Eq. B1):

```
−8πG/t = ⟨(2m² + t) P_J(1 + 2t/m²)⟩
```

where the average `⟨·⟩` is defined by the spectral integral in Eq. (7) of the main text. Because of this pole, the forward limit `t → 0` cannot be taken, and one must work with **fixed-t dispersion relations** at `t ≠ 0`, following Ref. [5] (Caron-Huot, Mazac, Rastelli, Simmons-Duffin).

### 2.2 The Functional Bootstrap Equation

Smearing both sides of the graviton-pole dispersion relation (Eq. B2) against a generic functional **α**, and using the spectrum ansatz (Eqs. 12 and 14), one arrives at the master equation (Eq. B3):

```
8πG · α·v_obj + λ_ϕ · α·v_norm + Σ_{i,J} λ_{m_i,J} · α·v_HE(m_i², J) + ⟨α·v_HE(m², J)⟩ = 0
```

The vectors entering this equation are defined component-by-component in **Eq. (B4)**:

- **v_obj = (2m₁^{1/2}, ..., 0, ...)**: the objective vector encoding the graviton residue;
- **v_HE(m², J) = (−G₁(m², J), ..., X_{1,0}(m², J), ...)**: the full dispersive+null constraint vector for a state of mass² `m²` and spin `J`, where G_k and X_{k,q} are the functions defined in Ref. [26];
- **v_norm = (−G₁(m²_ϕ, J_ϕ), ..., X_{1,0}(m²_ϕ, J_ϕ), ...)**: the same evaluated at the state being maximized.

### 2.3 The SDP Problem

One searches for **α** maximizing `α·v_obj` subject to (Eq. B5):

```
α·v_norm = 1
α·v_HE(m², J) ≥ 0   for all (m², J) in S
```

where `S` is the set of all isolated states below threshold plus all even-spin states above `M²`. This dual SDP problem, when solved, yields an upper bound on the normalized coupling `λ_ϕ`.

### 2.4 Numerical Parameters Cited in Appendix B

The paper lists the following default parameter values used throughout Section IV:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| J_max | 40 | Maximum spin in small-J region |
| J_huge | 5000 | Large-spin asymptotics cutoff |
| m_max | 10 | Maximum discrete mass above threshold |
| b_max | 80 | Maximum impact parameter |
| n_max | 13 | Maximum null constraint order |
| k_max | 10 | Maximum k in null constraints |

These are precisely reproduced in the code, as detailed below.

---

## 3. Code Structure: Global Setup

### 3.1 File and Parameter Input

```mathematica
Import["./SDPB2.m"]
xxx = ToExpression[ReadList["./range.txt"][[1]]]
mu1 = Rationalize[2/xxx]
mu2 = 1
mu3 = ToExpression[ReadList["./range.txt"][[3]]]
mu3 = Rationalize[mu3/xxx]
a1 = ToExpression[ReadList["./range.txt"][[2]]]
```

The code reads three parameters from an external file `range.txt`:
- **Entry 1 (`xxx`)**: a normalization/scaling factor.
- **Entry 2 (`a1`)**: used for bounding (±1 gives upper/lower bound search), though `a1` does not appear again explicitly in this file — it is presumably used when invoking SDPB externally.
- **Entry 3 (`mu3` before rescaling)**: the raw gap mass M.

After rescaling:
- `mu1 = 2/xxx` is the **mass² of the first (lightest) resonance**, chosen so that with `mu2 = 1` fixed (i.e., the second state has mass² = 1, providing the normalization scale), the ratio `mu1/mu2` is scanned by varying `xxx`.
- `mu3` is the **UV cutoff mass M²**, rescaled to units of `mu2`.

This parametrization directly implements the setup of Fig. 4 in the paper: the coupling `λ_{2,2} m₁⁶ / 8πG m₂²` is plotted as a function of `(m₁/m₂)²`, which is `mu1/mu2 = mu1` since `mu2 = 1`. Scanning `xxx` scans this ratio.

### 3.2 Global Dimensions

```mathematica
nulllist = {21, -1, -1, -1}
num1 = 27/2
```

- `nulllist[[1]] = 21`: the number of null constraints included. In Appendix B, `n_max = 13` for the null constraint order `n` (cf. paper's table of parameters), but this `nulllist` value counts the **total number of null constraint components** (i.e., the total size of the χ vector after summing over all valid `(n, k)` pairs for `n = 4, ..., n_max`).
- `num1 = 27/2`: determines the range of the smearing polynomials in the dispersive basis. The index `n` in the Cimp2 and Nlist functions runs from `3/2` to `num1` in steps of 1, giving `(27/2 - 3/2)/1 + 1 = 13` values — matching `n_max = 13` in the paper's table. This is the dimension of the smearing functional space in the dispersive sector.

```mathematica
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}]
list1 = Table[0, {i, 3/2, num1, 1}]
```

`list0` is a zero-padding array of length equal to the total number of null constraints (used to fill the null-constraint slots in v_HE for states that do not couple to certain blocks). `list1` is similarly a zero array over the dispersive sector.

---

## 4. Core Functional Kernels

### 4.1 The −2 Dispersive Kernel: `Cimp2list`

```mathematica
Cimp2list[n_, x_, J_] = -3 p^(1+n) x * (
    -2x Gamma[(1+n)/2] HypergeometricPFQRegularized[{-J,7+J,(1+n)/2},{4,(3+n)/2}, p²/x]
    + p² Gamma[(3+n)/2] HypergeometricPFQRegularized[{-J,7+J,(3+n)/2},{4,(5+n)/2}, p²/x]
) /. {p -> mu2}
```

This implements the **G_k(m², J) component** of the v_HE vector from Eq. (B4), specialized to the k = −2 improved subtraction (hence "Cimp2" = "C improved, 2 subtractions"). The arguments are:
- `n`: the moment index (runs from `3/2` to `num1` in the smearing);
- `x`: the mass² of the state (`m²`);
- `J`: the spin.

The function is evaluated at `p = mu2 = 1`, meaning the dispersion relations are fixed at `t = −p² = −1` (in units of the second-state mass). This is necessary because in the gravitational case one cannot take `t → 0`, and must work at fixed nonzero `p²`.

The Hypergeometric functions here arise from integrating the Gegenbauer polynomial `G_J^D(1 + 2t/m²)` against the subtraction kernel `1/[s(s−p²)]^{n/2+1}`, which in D=10 dimensions yields generalized hypergeometric functions of type `₃F₂`. The prefactor `−3 p^{1+n} x` accounts for the D=10 phase space factors and the −2 subtraction normalization. The term `(7+J)` in the hypergeometric parameters reflects the D=10 Gegenbauer algebra, where the relevant polynomial degree is `J + (D−4) = J + 6`, entering via the relation `{−J, (D−2)+J} = {−J, 7+J}` for D=10 (note D−2 = 8, but the hypergeometric parameters follow a specific convention tied to the Gegenbauer polynomial identity).

After evaluation at `p = mu2 = 1`, `Cimp2list[n, x, J]` is a rational function of `x` with `J`-dependent polynomial numerators, representing the contribution of a spin-J state of mass² `x` to the n-th moment of the dispersive sum rule.

### 4.2 The Null Constraint Kernel: `Nlist`

```mathematica
Nlist[n_, x_, J_] := x * { (polynomial in J and x for n=0), (polynomial for n=1), ... }[[n+1]]
```

`Nlist` (and its twin `Nlist2`) implement the **X_{k,q}(m², J) components** of the v_HE vector from Eq. (B4). These are the null constraint functionals derived from Eq. (8) of the main text,

```
χ_{n,k} = Res_{t=0}[ (1/s^{n−k−1}t^{k+1} − 1/t^{n−k−1}s^{k+1}) − (s → −s−t) ] P_ℓ(1 + 2t/s)
```

evaluated on the spectral density for a state of mass² `x` and spin `J`. The index `n` (running from 0 to `nulllist[[1]] = 21`) labels the null constraint, and the polynomial coefficients (which are rational functions of `x` with polynomial numerators in `J`) have been **pre-computed analytically** and hard-coded as a large explicit list. This is a substantial algebraic precomputation specific to D=10, encoding the action of 22 independent null constraint functionals (labeled `n = 0, ..., 21`) on a single-particle state.

The structure `x * {...}[[n+1]]` means: for each null constraint index `n`, the kernel is a degree-`(2n+2)` polynomial in `J` divided by a power of `x`, multiplied by an overall `x`. Many entries factor as `J(7+J) × polynomial(J)`, reflecting the fact that the null constraints vanish automatically for J=0 and J=−7 (which corresponds to J=0 in the shifted Gegenbauer convention for D=10), providing a useful cross-check.

`Nlist2` is a second version of the same function, used in the `PolyLargeJ` block, where the mass variable has been reparametrized for the large-J/large-mass regime. It contains identical polynomial coefficients to `Nlist`.

### 4.3 Large Impact Parameter Kernels: G, B, F, H

```mathematica
G[n_, D_, b_, mu_] := b^(num1+1) * mu^((1+n)/2) * (2^n Γ[(D−2)/2]Γ[(n+1)/2]) / Γ[(D−n−3)/2] / b^(n+1)
B[n_, D_, b_, mu_] := b^(num1+1) * mu^((1+n)/2) * (2^((D−3)/2) Γ[(D−2)/2]) / √π / b^((D−1)/2)
F[n_, D_, b_, mu_] := b^(num1+1) * mu^((1+n)/2) * (2^((D−9)/2) Γ[(D−2)/2]) / √π * (−27+12D−D²−8n) / b^((D−3)/2)
H[n_, D_, b_, mu_] := b^(num1+1) * mu^((1+n)/2) * Hypergeometric2F1[(1+n)/2,(D−2)/2,(3+n)/2, −mu*b²/4] / (n+1)
```

These four functions encode the **large impact parameter (large-b) asymptotics** of the dispersive functional basis, following Ref. [26] (Albert, Knop, Rastelli). In the impact parameter representation, the amplitude at large `b` is dominated by long-range gravitational exchange, and the functional constraints must encode positivity in this regime as well.

- **`G` and `B`** represent the leading power-law contributions from the spin-J expansion at large `b`, encoding the Bessel-function-type asymptotics of the Gegenbauer polynomials at large spin.
- **`F`** captures a subleading contribution proportional to `(−27 + 12D − D² − 8n) = (−27 + 120 − 100 − 8n) = (−7 − 8n)` for D=10; this is negative for all `n ≥ 0`, indicating a specific sign structure in the large-b falloff.
- **`H`** is the full Fourier-Bessel transform at fixed finite `b`, used for the discretized small-b region. It involves a `₂F₁` function whose argument `−mu*b²/4` grows with `b`.

The overall prefactor `b^(num1+1)` ensures the polynomial positivity structure required by the PMP format after factoring out the DampedRational prefactor.

All four are evaluated with `D = 10`, `mu = 1/mu2 = 1`, and various values of `b`.

---

## 5. Polynomial Matrix Blocks (Positivity Constraints)

The key mapping from the paper's Eq. (B5) to the code is through the `PositiveMatrixWithPrefactor` objects, each of which encodes `α · v_HE(m², J) ≥ 0` for a specific class of `(m², J)` values. Together they cover the full set `S` defined in Eq. (B5).

### 5.1 `PolySmall[J, x, y]`: Small Spin, Near-Threshold States

```mathematica
PolySmall[J_, x_, y_] := PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, y],
    {{ Flatten[Cancel[x^(Max[J, exp]) * {
        Table[Cimp2list[n, x, J], {n, 3/2, num1}],
        Table[Nlist[n, x, J], {n, 0, nulllist[[1]]}]
    }]] }}]
```

This block is used for **small even spins J = 0, 2, 4, ..., 40** and masses in the UV region (above threshold M, parametrized as `mu3 + x` with `x ≥ 0`). The positivity variable `y` is the mass variable above threshold. The `DampedRational[1, {}, 1/E, y]` prefactor provides an exponential damping `e^{−y}` which converts the semi-infinite mass integration into a polynomial-in-`y` positivity problem (the standard technique for SDPB).

The inner vector `{Cimp2list, Nlist}` is precisely the v_HE vector for a given (m², J): first the dispersive components (`Cimp2list` for `n = 3/2, ..., 27/2`), then the null constraint components (`Nlist` for `n = 0, ..., 21`). The prefactor `x^(Max[J, exp])` clears denominators in `x` (since `Cimp2list` and `Nlist` are rational in `x`), ensuring that the entries of the matrix are **polynomials** in `x` as required by the PMP format. The variable `exp = −Exponent[Nlist[21, x, 10], x, Min]` computes the most negative power of `x` appearing in the null constraints, so that `x^{exp}` is the minimal clearing prefactor.

### 5.2 `PolyMiddle[J, x, y]`: Middle Spin Range

```mathematica
PolyMiddle[J_, x_, y_] := PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, y],
    {{ Flatten[{
        Table[Cimp2list[n, x, J], {n, 3/2, num1}],
        Table[Nlist[n, x, J], {n, 0, nulllist[[1]]}]
    }] }}]
```

This is structurally identical to `PolySmall` but **without** the denominator-clearing prefactor `x^{Max[J,exp]}`. It is used for **middle spin values J = 42, 44, ..., 5000** (in various batches) and for the two isolated states below threshold:
- `PolyMiddle[0, mu2, x]`: spin-0 state at mass² = mu2 = 1 (the second resonance);
- `PolyMiddle[0, mu1, x]`: spin-0 state at mass² = mu1 (the first resonance).

For the isolated below-threshold states, the "mass variable" is fixed (a constant), and `x` here acts as the polynomial positivity variable of the `DampedRational` prefactor. The absence of the clearing prefactor is acceptable because for these applications the rational form of `Cimp2list`/`Nlist` in their mass argument is either trivially polynomial (for fixed mass) or the denominators are already handled by the mass parametrization.

Note that `PolyMiddle[0, mu1, x]` and `PolyMiddle[0, mu2, x]` encode the **below-threshold isolated state contributions** to the SDP, corresponding to the `λ_{m_i, J}` terms in Eq. (B3) for the two resonances in the two-state ansatz. The spin-2 and higher-spin couplings of these states are handled elsewhere (the spin-2 coupling of the second state is the one being maximized, entering v_norm and v_obj).

### 5.3 `PolyLargeJ[J, x, y]`: Large Spin Asymptotics

```mathematica
PolyLargeJ[J_, x_, y_] := PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, y],
    {{ Flatten[Simplify[{
        Table[0, {n, 3/2, num1}],
        Table[(1*mu3 + 10*mu3*x)^exp * Nlist2[n, (1*mu3 + 10*mu3*x)/(1+x), J],
              {n, 0, nulllist[[1]]}]
    }]] }}]
```

This block handles the **large-J continuous limit** (Ref. [26]). At large spin, the discrete sum over J can be replaced by a continuous integral, and the SDP constraint becomes `α · v_HE(m², J) ≥ 0` for continuous `J`. This is encoded as a polynomial-positivity constraint in the variable `x` (which here parametrizes the spin J, not the mass).

Key observations:
1. **The dispersive (Cimp2) components are set to zero** — at very large J the dispersive contribution is suppressed and the null constraints dominate.
2. **Mass reparametrization**: `m² = (mu3 + 10·mu3·x)/(1+x)` maps `x ∈ [0, ∞)` to `m² ∈ [mu3, 10·mu3]`, covering the UV region above threshold M up to `10M`. This compactification ensures polynomial behavior in the positivity variable `x`.
3. **`Nlist2`** (identical in content to `Nlist`) is used here rather than `Nlist`, presumably for clarity or to avoid symbol conflicts in parallel evaluation.
4. The prefactor `(mu3 + 10*mu3*x)^exp` again clears the denominator in the mass variable.

The usage in `TSDP` covers `J` ranges of {10000–50000 step 5000} and {J+5000} with varying `b`, matching the paper's `J_huge = 5000` parameter.

### 5.4 `NewPoly[b]`: Large Impact Parameter (Continuous b)

```mathematica
NewPoly[b_] := PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, x],
    {
      { G+B, F },
      { F,   G−B }
    }]
```

where each matrix entry is `Flatten[{Table[G[i,10,b,1/mu2], ...], list0}]` etc. This block encodes the **large-impact-parameter positivity constraint** as a **2×2 positive semidefinite matrix**. The structure

```
M = [ G(b) + B(b),   F(b) ]
    [ F(b),          G(b) − B(b) ]
```

must be PSD, which requires both diagonal entries ≥ 0 and `det(M) = (G+B)(G−B) − F² = G² − B² − F² ≥ 0`. This 2×2 structure encodes the combination of **even and odd-in-b contributions** to the large-b asymptotics of the functional constraints. The off-diagonal `F` term mixes them, implementing the constraint structure of Ref. [26] adapted to the gravitational case.

The `list0` zero-padding fills the null-constraint slots (large-b constraints live entirely in the dispersive sector). The variable `x` is the positivity variable, and `b` is treated as a **fixed continuous parameter** at the value `b = 80` (the upper boundary `NewPoly[x−80]` effectively enforces `b ≤ 80`, i.e., b_max = 80 from the paper's table).

### 5.5 `Polydis[b, x]`: Discretized Small Impact Parameter

```mathematica
Polydis[b_, x_] := PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, x],
    {{ Flatten[{Table[H[n, 10, b, 1/mu2], {n, 3/2, num1}], list0}] }}]
```

This covers **small and medium b values discretely** (b from 1/32 to 80 in steps of 1/32), using the full hypergeometric kernel `H` which is valid at any b (not just asymptotically). The single-row matrix (1×1 PSD = non-negative) enforces `α · v_HE(b) ≥ 0` at each discrete value of `b`. The range `b ∈ [1/32, 80]` with step `1/32` gives 2560 discrete b-points, providing dense coverage of the impact parameter domain.

---

## 6. The Norm and Objective Vectors

```mathematica
norm = 1 * N[Flatten[{
    Table[Cimp2list[i, mu1, 2], {i, 3/2, num1}],
    Table[Nlist[n, mu1, 2], {n, 0, nulllist[[1]]}]
}], 600]

obj = −1 * N[Flatten[{
    Table[1/(i−1) * mu2^((1−i)/2), {i, 3/2, num1}],
    list0
}], 600]
```

These correspond directly to **v_norm** and **v_obj** from Eq. (B4):

**`norm` ↔ v_norm**: This evaluates the full v_HE vector at the state being maximized: mass² = `mu1` (the first resonance mass²), spin J = 2 (for the two-state problem bounding λ_{2,2}). For the three-state problem bounding λ_{3,4}, one would use `mu3` and J=4 instead. The normalization condition `α · v_norm = 1` (Eq. B5 first line) is implemented by passing `norm` to the SDPB problem.

**`obj` ↔ v_obj**: The objective vector has components `1/(i−1) · mu2^{(1−i)/2}` for `i = 3/2, 5/2, ..., 27/2`. Since `mu2 = 1`, this simplifies to `2/(2i−2) = 1/(i−1)`. The paper's Eq. (B4) gives `v_obj = (2m₁^{1/2}, 0, 0, ...)` — the `2m₁^{1/2}` coefficient multiplies the leading smearing basis function. In the Laguerre/polynomial smearing basis parametrized by `i`, the expansion `2√m₁ = Σ_i c_i · φ_i(m₁)` produces these rational coefficients `1/(i−1)`, which arise from the Laguerre polynomial expansion of the function `2√m at m = mu1`. The null-constraint slots are zero (`list0`), since the graviton contributes only to the dispersive sector. The overall sign `−1` is because SDPB minimizes the objective; we want to **maximize** `α · v_obj`, so we **minimize** `−α · v_obj`.

The 600-digit floating-point precision (`N[..., 600]`) is required because SDPB works in arbitrary precision arithmetic, and the coefficients of `Nlist` involve very large integers (as seen in the code); loss of precision here would corrupt the SDP solution.

---

## 7. Assembly: The `TSDP` Function

```mathematica
TSDP[datfile_] := Module[{
    pols = Flatten[{
        (* 1 *) PolySmall for J=0..40 above threshold
        (* 2 *) PolyMiddle for J=0 at mu2 (second resonance)
        (* 3 *) PolyMiddle for J=0 at mu1 (first resonance)
        (* 4–6 *) PolyMiddle for J=42..5000 above threshold
        (* 7–9 *) PolyLargeJ for large J
        (* 10 *) NewPoly at b=80 (boundary)
        (* 11 *) Polydis for b=1/32..80
        (* 12 *) More PolyMiddle J=5500..10000
        (* 13 *) More PolyLargeJ J=55000..100000
    }, 1],
    norm = ...,
    obj = ...
},
WritePmpJson[datfile, SDP[obj, norm, pols]]]
```

The `pols` list is the complete set of positivity constraints encoding `α · v_HE(m², J) ≥ 0` for all `(m², J) ∈ S`. The assembly is organized into 13 batches, which together cover the following domains:

| Batch | Function | J range | Mass range | Purpose |
|-------|----------|---------|------------|---------|
| 1 | `PolySmall` | 0, 2, ..., 40 | > mu3 (UV) | Small even spin above M |
| 2 | `PolyMiddle` | 0 | mu2 = 1 | Second resonance, spin-0 |
| 3 | `PolyMiddle` | 0 | mu1 | First resonance, spin-0 |
| 4 | `PolyMiddle` | 42–100 step 2 | 1/j, j ∈ 1/1000mu3..1/mu3 | Medium J above M |
| 5 | `PolyMiddle` | 100–900 step 50 | same | Higher medium J above M |
| 6 | `PolyMiddle` | 1000–5000 step 500 | same | High J above M |
| 7 | `PolyLargeJ` | 10000–50000 step 5000 | compactified | Very large J continuous |
| 8 | `PolyLargeJ` | J+5000, b=10..300 | compactified | Large J, large b |
| 9 | `PolyLargeJ` | J+5000, b=0..10 | compactified | Large J, small b |
| 10 | `NewPoly` | — | — | b = 80 boundary |
| 11 | `Polydis` | — | — | b = 1/32..80, discrete |
| 12 | `PolyMiddle` | 5500–10000 step 500 | same as 4–6 | Extended J coverage |
| 13 | `PolyLargeJ` | 55000–100000 step 5000 | compactified | Extended large J |

The parallel table evaluation (`ParallelTable` with `LaunchKernels[]`) distributes the heavy polynomial arithmetic across available CPU cores, since computing the full `PolySmall`/`PolyMiddle` entries at 600-digit precision is computationally expensive.

The final call `WritePmpJson["out.json", SDP[obj, norm, pols]]` formats the entire SDP as a JSON file in the SDPB PMP format, ready for submission to the solver.

---

## 8. Correspondence Summary Table

| Code Element | Paper Reference | Physical Meaning |
|---|---|---|
| `mu1 = 2/xxx` | m₁² (first resonance mass²) | Scanned parameter in Figs. 4, 5 |
| `mu2 = 1` | m₂² (second resonance mass²) | Normalization scale |
| `mu3` | M² (UV gap/cutoff) | Threshold above which states are continuous |
| `nulllist = {21,...}` | n_max, k_max in Eq. (8) | 22 null constraint components |
| `num1 = 27/2` | n_max = 13 smearing moments | Dimension of dispersive functional |
| `Cimp2list[n,x,J]` | G_k(m²,J), Eq. (B4) | Dispersive sum rule kernel, −2 subtraction |
| `Nlist[n,x,J]` | X_{k,q}(m²,J), Eq. (B4) | Null constraint kernel |
| `G, B, F, H` | Large-b asymptotics, Ref.[26] | Impact parameter positivity |
| `PolySmall/Middle` | v_HE(m²,J) ≥ 0, Eq. (B5) | Positivity for discrete (m²,J) |
| `PolyLargeJ` | v_HE(m²,J) ≥ 0 continuous J | Large-J continuum positivity |
| `NewPoly` | Large-b 2×2 PSD constraint | Impact parameter boundary |
| `Polydis` | Small-b positivity | Finite-b sum rule constraint |
| `norm` | v_norm, Eq. (B4) | Normalization of coupling being bounded |
| `obj` | v_obj, Eq. (B4) | Graviton contribution / objective |
| `SDP[obj, norm, pols]` | Eqs. (B3)–(B5) | Full dual SDP problem |
| `WritePmpJson` | SDPB input format | Output for SDPB solver |

---

## 9. Key Technical Observations

### 9.1 Rationalization and Precision

The use of `Rationalize[mu1]` and `Rationalize[mu3]` ensures that the mass² values entering `Cimp2list` and `Nlist` are **exact rational numbers**, not floating-point approximations. This is crucial because `Nlist` contains rational functions of the mass with large integer coefficients (up to ~10^29 in the numerators), and any floating-point error in the mass would propagate through all those terms. The 600-digit precision in `N[..., 600]` provides the final numerical output at sufficient precision for SDPB.

### 9.2 The Exponent Clearing Strategy

```mathematica
exp = -Exponent[Nlist[nulllist[[1]], x, 10], x, Min]
```

This computes the most negative power of x appearing in the null constraint list (at the highest n and J=10 as a representative value). The result `exp` is then used in `PolySmall` to multiply the entire vector by `x^{Max[J,exp]}`, clearing all negative powers of x. This converts rational functions of x into polynomials, which is the input format required by SDPB's PMP. The `Max[J, exp]` handles both the case where the null constraints have the most severe denominators (taking `exp`) and the case where the dispersive term with high J has `J` powers of x in the numerator (taking J).

### 9.3 The Two-State vs. Three-State Structure

The current code targets the **two-state problem** (Fig. 4 in the paper): two resonances below threshold M, with the spin-2 coupling λ_{2,2} at mass mu2 being maximized. This is encoded by:
- `PolyMiddle[0, mu2, x]` and `PolyMiddle[0, mu1, x]` as the only two below-threshold isolated states (both with spin J=0; higher spins of these states are absent in the two-state ansatz by construction);
- `norm` evaluated at (mu1, J=2): the state whose coupling is normalized.

For the **three-state problem** (Fig. 5), one would additionally include `PolyMiddle[J, mu3_below, x]` blocks for the third state, and change `norm` to evaluate at (mu3, J=4). The code structure is clearly designed to accommodate this extension.

### 9.4 The Graviton Constraint

The graviton pole does not appear explicitly as a separate polynomial block, because it is absorbed into the objective and norm vectors through Eqs. (B1)–(B4). Specifically, the graviton's contribution to the dispersion relation is encoded in the linear equality `α · v_norm = 1` (norm vector) and the gravitational coupling `8πG` appears through the ratio in the objective normalization `λ_{2,2} m₁⁶ / 8πG m₂²` visible in Fig. 4's y-axis label. The `obj` vector's `1/(i−1) mu2^{(1−i)/2}` entries represent the expansion of `2/√m₁` (the graviton pole residue's m₁ dependence) in the smearing basis.

---

## 10. Workflow Summary

The complete computational pipeline is:

1. **`range.txt`**: Provides `(xxx, a1, mu3_raw)`, specifying the mass ratio `(m₁/m₂)²` point being computed and the cutoff M.
2. **`generate.m`**: Reads the parameters, constructs all polynomial positivity blocks (13 batches), assembles `norm` and `obj`, writes `out.json`.
3. **SDPB** (external): Reads `out.json`, solves the dual SDP, returns the optimal `α · v_obj` as the upper bound on the normalized coupling.
4. **Scanning**: Steps 1–3 are repeated for a grid of `xxx` values, producing the curves shown in Figs. 4 and 5 of the paper.

Each `xxx` value corresponds to one horizontal data point in those figures. The parameter `a1 = ±1` controls whether an upper or lower bound is sought (though in the figures only upper bounds are shown).

---

*Report prepared by analysis of arXiv:2510.07991v1 and the accompanying Mathematica source file `generate.m`.*
