(* Selecting Basis *)

Clear[basisFunctions];

(* For D = odd and D >= 7, we can use a simple set of basis functions that capture the essential behavior. *)

basisFunctions[D_, maxPow_] := Module[{list, k},
  (* include p^2, p^3, and p^n in 2 ... N_f *)
  list = Table[p^n , {n, 2, N_f}];
  list
];

(* Precompute Analytic Kernel *)

(* We can compute the kernel analytically for the chosen basis functions. *)

T1[n_, J_, m_]:= Integrate[ p^n * (2m^2 - p^2) * PJ[J, 1- 2*p^2/m^2] / ( m^2*(m^2 - p^2)^2 ) ,{p,0,1}];

(* Calculation of derivative of Ge *)

(* EFT action integral *)

gamma[n_]:= Integrate[p^n / p^2, {p,0,1}];

alpha[n_]:= Integrate[p^n, {p,0,1}];

beta[n_]:= Integrate[p^n * p^2, {p,0,1}];

(* ===================================================================
   Null Constraint Integrals
   Source: eq:nullconstraints in sec_localized_b.tex
   ===================================================================

   The null constraint kernel X_k[m^2, J, u] is defined by:

     X_k[m^2,J,u] =
       (2m^2+u) * m^2 * P_J(1+2u/m^2)
       ------------------------------------
         (u * m^2 * (m^2+u))^{1+k/2}

       -  Res_{u'=0} [
            (2m^2+u')(m^2-u')(m^2+2u')         m^2 * P_J(1+2u'/m^2)
            --------------------------------- * ----------------------
            m^2(u-u')u'(m^2-u)(m^2+u')(m^2+u+u')  (u'*m^2*(m^2+u'))^{k/2}
          ]

   where k = 4, 6, ..., and P_J is the Gegenbauer polynomial normalized
   to P_J(1) = 1, i.e. P_J(x) = 2F1(-J, J+D-3; (D-2)/2; (1-x)/2).

   Substituting u = -p^2 gives the kernel integrated against h_k(p) = sum_i b_{k,i} p^i.
   The integral NullConstraintIntegral[k, i, J, m] = integral_0^1 dp p^i * X_k(m^2,J,-p^2)
   is "KernelX[k][i][J](m)" in the pseudocode (Table B).

   Key properties:
     - <X_k> = 0 for any UV-complete spectrum (null constraint)
     - The residue subtracts all EFT contributions, leaving a purely UV quantity
     - The pole order at u'=0 is 1 + k/2 (order 3 for k=4, order 4 for k=6)
=================================================================== *)

(* Normalized Gegenbauer polynomial P_J(x) with P_J(1) = 1.
   Paper convention: P_J(x) = 2F1(-J, J+D-3; (D-2)/2; (1-x)/2)
   In Mathematica: GegenbauerP[J, lambda, x] with lambda = (D-3)/2,
   divided by its value at x=1 to enforce the normalization P_J(1) = 1.
   Check: P_J'(1) = J*(J+D-3)/(D-2)  [quadratic Casimir / (D-2)].
   Note: D is a global parameter set before calling these functions. *)

PJ[J_, x_] :=
  GegenbauerP[J, (D - 3)/2, x] / GegenbauerP[J, (D - 3)/2, 1];

(* First (direct) term of X_k at general u.
   Combines the two factors from the paper into a single rational expression:
     (2m^2+u) / (u*m^2*(m^2+u))  *  m^2*P_J(1+2u/m^2) / (u*m^2*(m^2+u))^{k/2}
   = (2m^2+u) * P_J(1+2u/m^2) / (u*m^2*(m^2+u))^{1+k/2}               *)

XKernelFirstTerm[k_, m_, J_, u_] := (2*m^2 + u) * (m^2) * PJ[J, 1 + 2*u/m^2] / (u * m^2 * (m^2 + u))^(1 + k/2);

(* Residue integrand: the expression inside Res_{u'=0} in eq:nullconstraints.
   The pole at u'=0 has order 1 + k/2, arising from:
     - one explicit factor of u' in the denominator
     - k/2 factors from (u'*m^2*(m^2+u'))^{k/2}
   Mathematica's Residue[] computes the Laurent coefficient of u'^{-1}.  *)

XKernelResidueIntegrand[k_, m_, J_, u_, up_] :=
  (2*m^2 + up) * (m^2 - up) * (m^2 + 2*up) /
  (m^2 * (u - up) * up * (m^2 - u) * (m^2 + up) * (m^2 + u + up)) *
  PJ[J, 1 + 2*up/m^2] / (up * m^2 * (m^2 + up))^(k/2);

(* Full null constraint kernel X_k[m^2, J, u]:
   first term minus the residue subtraction at u'=0.
   When evaluated at u = -p^2 (spacelike momentum transfer), this is
   the kernel that multiplies h_k(p) in the positivity constraint. *)

XKernel[k_, m_, J_, u_] :=
  XKernelFirstTerm[k, m, J, u] -
  Residue[XKernelResidueIntegrand[k, m, J, u, up], {up, 0}];

(* Integral of p^i * X_k[m^2, J, -p^2] over p in [0, 1].
   This is the entry KernelX[k][i][J](m) from the precomputation table.
   It is the contribution of the basis function h_k(p) = p^i to the
   positivity constraint F(m, J) >= 0 at mass m and spin J.

   Usage:  NullConstraintIntegral[4, 0, J, m]  -- k=4, i=0 (constant h_4)
           NullConstraintIntegral[4, 1, J, m]  -- k=4, i=1 (linear h_4)
           NullConstraintIntegral[6, 0, J, m]  -- k=6, i=0  etc.

   Self-check (A-6 from research.md): for a test spectral density rho_J(m),
   summing over (m, J) with the UV measure should give zero, since <X_k> = 0. *)

NullConstraintIntegral[k_, i_, J_, m_] :=
  Integrate[p^i * XKernel[k, m, J, -p^2], {p, 0, 1},
    Assumptions -> {m > 1, k > 0, i >= 0, Element[J, Integers], J >= 0}];


(* I DON'T THINK WE NEED IMPACT PARAMETER SPACE INTEGRALS *)


(* ===================================================================
   Impact Parameter Space Integrals
   Source: app_flat_space_numerics.tex (eq:the1f2s) and
           sec_localized_b.tex (flat fourier / eq:impactparametertransform)
   ===================================================================

   In the m -> infinity limit with fixed impact parameter b = 2J/m,
   the positivity constraint on F(m, J) becomes:

     Gamma((D-2)/2) * integral_0^1 dp f(p) * J_{(D-4)/2}(b*p) / (b*p/2)^{(D-4)/2} >= 0

   for all b >= 0.  The null-constraint functions h_k(p) are subleading
   in this m -> infinity limit and do not enter (Block-2 in the pseudocode).

   For a pure power f(p) = p^n the integral evaluates analytically to:

     Gamma((D-2)/2) * integral_0^1 dp p^n * J_{(D-4)/2}(b*p) / (b*p/2)^{(D-4)/2}
       = _1F_2( (n+1)/2 ; (D-2)/2, (n+3)/2 ; -b^2/4 ) / (n+1)

   This is "ImpactExact[n](b)" in the pseudocode (Table D).

   Self-check (A-3 from research.md):
     ImpactIntegral[n, 0] = 1/(n+1) = alpha[n]
   because _1F_2(...; 0) = 1 by the series definition.
=================================================================== *)

(* Exact analytic impact parameter integral for basis function p^n.
   Valid for all b >= 0, any integer or half-integer n > 1.
   At b = 0: reduces to 1/(n+1) = alpha[n] (the g_2 action integral).  *)

ImpactIntegral[n_, b_] :=
  HypergeometricPFQ[{(n + 1)/2}, {(D - 2)/2, (n + 3)/2}, -b^2/4] / (n + 1);

(* Verification of self-check A-3: ImpactIntegral[n, 0] should equal alpha[n] = 1/(n+1).
   This holds because HypergeometricPFQ[{a},{b,c},0] = 1 for any a,b,c. *)

(* alpha[n] = 1/(n+1), so ImpactIntegral[n, 0] = 1/(n+1) = alpha[n]. Confirmed. *)


(* Large-b asymptotic expansion of ImpactIntegral[n, b] (eq:the1f2s):
   The _1F_2 function has the asymptotic form:

     _1F_2((n+1)/2; (D-2)/2, (n+3)/2; -b^2/4) / (n+1)
       ~  A_n / b^{n+1}
        + B_n * cos(b - pi*(D-1)/4) / b^{(D-1)/2}
        + C_n * sin(b - pi*(D-1)/4) / b^{(D-1)/2 + 1}
        + ...

   where the non-oscillatory coefficient is:
     A_n = 2^n * Gamma((D-2)/2) * Gamma((n+1)/2) / Gamma((D-n-3)/2)

   and the leading oscillatory coefficient is n-independent:
     B_n = 2^{(D-3)/2} * Gamma((D-2)/2) / Sqrt[Pi]

   A_n governs whether the impact-parameter positivity can be achieved:
   if n > (D-3)/2, the non-oscillatory term is subleading relative to
   the oscillatory term at large b, making positivity hard to enforce
   (this is why D=5 requires the difference basis p^k - p^2).           *)

(* Non-oscillatory coefficient of 1/b^{n+1} term *)
ImpactCoeffA[n_] :=
  2^n * Gamma[(D - 2)/2] * Gamma[(n + 1)/2] / Gamma[(D - n - 3)/2];

(* Leading oscillatory coefficient of cos(b - pi*(D-1)/4) / b^{(D-1)/2} term.
   This coefficient is the same for all n at leading order in 1/b.      *)
ImpactCoeffB[n_] :=
  2^((D - 3)/2) * Gamma[(D - 2)/2] / Sqrt[Pi];

(* Truncated large-b asymptotic form, keeping only the two leading terms.
   Used to build the Block-3 PSD matrix (eq:strongercondition):
     [[A(b)+B(b), C(b)], [C(b), A(b)-B(b)]] >= 0  for b >= B_large
   where A(b) = sum_n a_n * ImpactCoeffA[n]/b^{n+1},
         B(b) = sum_n a_n * ImpactCoeffB[n] * cos-coefficient / b^{(D-1)/2},
         C(b) = sum_n a_n * ImpactCoeffB[n] * sin-coefficient / b^{(D-1)/2+1}. *)

ImpactAsymp[n_, b_] :=
  ImpactCoeffA[n] / b^(n + 1) +
  ImpactCoeffB[n] * Cos[b - Pi*(D - 1)/4] / b^((D - 1)/2);