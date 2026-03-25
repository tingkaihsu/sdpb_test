ClearAll["Global`*"]


(* ===================================================================
   Gravitational EFT Bootstrap: g2-g3 Bounds via S-Matrix Positivity
   Reference: Caron-Huot, Mazac, Rastelli, Simmons-Duffin
   Units: M = 1, 8*Pi*G = 1 throughout (paper convention, App. A)
   =================================================================== *)

(* ===================================================================
   GLOBAL PARAMETERS (set before calling any function)
   =================================================================== *)

d  = 5;     (* spacetime dimension; change to 5,6,... as needed *)
Nf = 7;     (* number of basis functions for f(p);
               d=5: p^3-p^2,...,p^{Nf+2}-p^2 = {p^3-p^2,...,p^9-p^2} *)

(* ===================================================================
   SECTION 1: Gegenbauer Polynomial, Normalized to P_J(1) = 1

   Paper convention (B2 improved flat, sec_localized_b.tex):
     P_J(x) = 2F1(-J, J+d-3; (d-2)/2; (1-x)/2)
   satisfies P_J(1) = 1 for all J.

   In Mathematica, GegenbauerP[J, lambda, x] with lambda = (d-3)/2
   is the standard Gegenbauer C^lambda_J(x).  It equals the above
   up to the overall value C^lambda_J(1), so we divide by it.

   Derivative identity:
     P_J'(1) = J*(J+d-3)/(d-2)
   This is the ONLY value of the derivative that ever appears in
   the kernel C2imp (in the T3 term below).  We do NOT need the
   full derivative function.
   =================================================================== *)

PJ[J_Integer, x_?NumericQ] :=
  GegenbauerP[J, (d-3)/2, x] / GegenbauerP[J, (d-3)/2, 1];

(* PJprime and m^2 = 1 *)
PJprime1[J_] := J*(J + d - 3) / (d - 2);

(* ===================================================================
   THIS IS CURRENTLY WRONG USING WRONG BASE FUNCTIONS IN d = 5
   BUT FIX IN SECTION 4.
   SECTION 2: Improved Kernel C2imp and Its Three Sub-Integrals

   From eq. (B2 improved flat) in sec_localized_b.tex, with u = -p^2:

     C2imp[m^2, J, u] =
         (2m^2+u)*P_J(1+2u/m^2)           u^2   (4m^2+3u)*P_J(1)   4u*P_J'(1)
         -------------------------    -   ---- * [---------------  + ----------]
            m^2*(m^2+u)^2                 m^6       (m^2+u)^2          m^4-u^2

   Split into three terms:

     T1: the first fraction (J-dependent via P_J)
     T2: the (4m^2+3u)*P_J(1) piece (J-INDEPENDENT since P_J(1)=1)
     T3: the 4u*P_J'(1)/(m^4-u^2) piece (J-dependent via P_J'(1))

   Substituting u = -p^2 in each:

     T1 integrand: p^n*(2m^2-p^2)*P_J(1-2p^2/m^2) / [m^2*(m^2-p^2)^2]
     T2 integrand: -p^n * p^4*(4m^2-3p^2) / [m^6*(m^2-p^2)^2]
     T3 integrand: +p^n * 4*p^6*P_J'(1) / [m^6*(m^4-p^4)]

   Sign derivation for T3:
     The term is -(u^2/m^6)*4u*P_J'(1)/(m^4-u^2)
     Under u=-p^2: u^2=p^4, 4u=-4p^2, m^4-u^2=m^4-p^4
     So: -(p^4/m^6)*(-4p^2)*P_J'(1)/(m^4-p^4)
       = +4*p^6*P_J'(1)/(m^6*(m^4-p^4))  [POSITIVE sign]
   =================================================================== *)

(* T1[n_, J_, m_] := Integrate[
  p^n * (2*m^2 - p^2) * PJ[J, 1 - 2*p^2/m^2] /
  (m^2 * (m^2 - p^2)^2),
  {p, 0, 1},
  Assumptions -> {m > 1, n > 1}
]; *)

T1[n_, J_, m_?NumericQ] := NIntegrate[
  p^n * (2*m^2 - p^2) * PJ[J, 1 - 2*p^2/m^2] /
  (m^2 * (m^2 - p^2)^2),
  {p, 0, 1},
  WorkingPrecision -> 50
];

T2[n_, m_?NumericQ] := NIntegrate[
  -p^n * p^4 * (4*m^2 - 3*p^2) /
  (m^6 * (m^2 - p^2)^2),
  {p, 0, 1},
  WorkingPrecision -> 50
];

(* T3 factored: pull PJprime1[J] outside the integral since it is
   a number at fixed J.  The integral itself is J-independent. *)
T3[n_, J_, m_?NumericQ] := PJprime1[J] * NIntegrate[
  4 * p^n * p^6 /
  (m^6 * (m^4 - p^4)),
  {p, 0, 1},
  WorkingPrecision -> 50
];

(* Full kernel integral for pure power p^n: *)
KernelC[n_, J_, m_?NumericQ] := T1[n, J, m] + T2[n, m] + T3[n, J, m];

(* ===================================================================
   SECTION 3: EFT Action Integrals

   From eq. (fulllinearprogramagain):
     I[f] = Int_0^1 dp f(p) * [1/p^2 + 2*g2 + g3*p^2]
           = gamma[f] + 2*g2*alpha[f] + g3*beta[f]

   For pure power p^n (analytical, no integration needed):
     gamma[n] = 1/(n-1)   [requires n > 1]
     alpha[n] = 1/(n+1)
     beta[n]  = 1/(n+3)

   Self-check: these equal the naive integrals directly, consistent
   with alpha[n] = ImpactIntegral[n, 0] (A-3 check in research-1.md).
   =================================================================== *)

gammaF[n_] := 1/(n - 1);
alphaF[n_] := 1/(n + 1);
betaF[n_]  := 1/(n + 3);

(* ===================================================================
   SECTION 4: d=5 Difference Basis p^k - p^2

   From Table 1 (tab:basisfunctionsforf) in sec_localized_b.tex:
   For d=5, basis = {p^3-p^2, p^4-p^2, ..., p^{Nf+2}-p^2}.

   Reason: in d=5, (d-3)/2 = 1 < 2, so pure powers p^n (n>=2)
   have their large-b Fourier transform dominated by the oscillatory
   term ~cos(b-pi*(d-1)/4)/b^{(d-1)/2} = cos(b-pi)/b^2, which
   cannot be made globally positive.
   The B_n coefficient (see eq:the1f2s) is n-independent:
     B_n = 2^{(d-3)/2}*Gamma((d-2)/2)/sqrt(pi)
   so the difference (p^k - p^2) exactly cancels the oscillatory
   B-coefficient, leaving only the non-oscillatory A_k/b^{k+1} term.

   For Nf=7: basisPowers = {3,4,5,6,7,8,9}
   Each basis function: basisFn[[i]](p) = p^{basisPowers[[i]]} - p^2

   Action integrals for the difference basis:
     gamma_d[k] = gamma[k] - gamma[2] = 1/(k-1) - 1
     alpha_d[k] = alpha[k] - alpha[2] = 1/(k+1) - 1/3
     beta_d[k]  = beta[k]  - beta[2]  = 1/(k+3) - 1/5
   =================================================================== *)

basisPowers = Range[3, Nf + 2];   (* {3,4,5,6,7,8,9} for Nf=7 *)

gamma_d[k_] := gammaF[k] - gammaF[2];   (* 1/(k-1) - 1 *)
alpha_d[k_] := alphaF[k] - alphaF[2];   (* 1/(k+1) - 1/3 *)
beta_d[k_]  := betaF[k]  - betaF[2];    (* 1/(k+3) - 1/5 *)

(* C_Improved *)
(* Kernel integral for the difference basis function (p^k - p^2): *)
KernelC_d[k_, J_, m_] := KernelC[k, J, m] - KernelC[2, J, m];

(* ===================================================================
   SECTION 5: Null Constraint Kernel X_k and Its Integrals

   From eq. (nullconstraints) in sec_localized_b.tex:

     X_k[m^2, J, u] =
         (2m^2+u) * m^2 * P_J(1+2u/m^2)
         ----------------------------------- (1)
           (u * m^2 * (m^2+u))^{1+k/2}

       - Res_{u'=0} [
           (2m^2+u')(m^2-u')(m^2+2u')        m^2 * P_J(1+2u'/m^2)
           --------------------------------- * -------------------------
           m^2*(u-u')*u'*(m^2-u)*(m^2+u')    (u'*m^2*(m^2+u'))^{k/2}
           *(m^2+u+u')
         ]

   (1) Note: the paper writes the first term as:
     [(2m^2+u)/(u*m^2*(m^2+u))] * [m^2*P_J(...)/(u*m^2*(m^2+u))^{k/2}]
   which combines to (2m^2+u)*m^2*P_J(...) / (u*m^2*(m^2+u))^{1+k/2}.
   We keep m^2 explicit in the numerator to match the paper's formula exactly.

   The residue at u'=0 has pole order 1+k/2 (order 3 for k=4, 4 for k=6).
   Mathematica's Residue[] extracts the coefficient of u'^{-1} in the
   Laurent expansion, which is exactly what we need.

   h_k basis: integer powers p^i, i = 0,1,...,N_h[k]-1
   For Figure 1: N_h[4]=6 (i=0..5), N_h[6]=4 (i=0..3)
   ===================================================================*)

(* what is NhByK??? *)
(* This is a Mathematica Association (a key-value dictionary). The syntax <| key -> value, … |> constructs it. *)
NhByK = <|4 -> 6, 6 -> 4|>;

XKernelFirstTerm[k_, m_, J_, u_] :=
  (2*m^2 + u) * m^2 * PJ[J, 1 + 2*u/m^2] /
  (u * m^2 * (m^2 + u))^(1 + k/2);

XKernelResidueIntegrand[k_, m_, J_, u_, up_] :=
  ((2*m^2 + up) * (m^2 - up) * (m^2 + 2*up)) /
  (m^2 * (u - up) * up * (m^2 - u) * (m^2 + up) * (m^2 + u + up)) *
  PJ[J, 1 + 2*up/m^2] / (up * m^2 * (m^2 + up))^(k/2);

XKernel[k_, m_, J_, u_] :=
  XKernelFirstTerm[k, m, J, u] -
  Residue[XKernelResidueIntegrand[k, m, J, u, up], {up, 0}];

(* KernelX[k, n, J, m]: integral of p^n * X_k[m^2,J,-p^2] over p in [0,1].
   This is "KernelX[k][n][J](m)" in the pseudocode Table B. *)
KernelX[k_, n_, J_, m_] :=
  Integrate[
    p^n * XKernel[k, m, J, -p^2],
    {p, 0, 1},
    Assumptions -> {m > 1, k > 0, n >= 0, Element[J, Integers], J >= 0}
  ];

(* ===================================================================
   SECTION 6: Impact Parameter Space Integrals
   
   Not sure what this is for?

   From sec:impactparamineqs in app_flat_space_numerics.tex:
   The m->infinity limit of the kernel positivity condition at fixed
   impact parameter b = 2J/m gives:

     Gamma((d-2)/2) * Int_0^1 dp p^n * J_{(d-4)/2}(b*p) / (b*p/2)^{(d-4)/2}
       = _1F_2( (n+1)/2 ; (d-2)/2, (n+3)/2 ; -b^2/4 ) / (n+1)

   Self-check A-3: at b=0, _1F_2(...;0)=1, so ImpactIntegral[n,0]
                   = 1/(n+1) = alpha[n].  Confirmed trivially.

   For the difference basis p^k - p^2:
   ImpactIntegral_d[k,b] = ImpactIntegral[k,b] - ImpactIntegral[2,b]
   =================================================================== *)

ImpactIntegral[n_, b_] :=
  HypergeometricPFQ[{(n+1)/2}, {(d-2)/2, (n+3)/2}, -b^2/4] / (n+1);

ImpactIntegral_d[n_, b_] :=
  ImpactIntegral[n, b] - ImpactIntegral[2, b];

(* Large-b asymptotic coefficients (eq:the1f2s):
     ImpactIntegral[n,b] ~ A_n/b^{n+1} + B_n*cos(b-pi*(d-1)/4)/b^{(d-1)/2} + ...
   A_n (non-oscillatory):  *)
ImpactCoeffA[n_] :=
  2^n * Gamma[(d-2)/2] * Gamma[(n+1)/2] / Gamma[(d-n-3)/2];

(* B_n (leading oscillatory): same for ALL n at leading order *)
ImpactCoeffB[n_] :=
  2^((d-3)/2) * Gamma[(d-2)/2] / Sqrt[Pi];

(* For the difference basis in d=5:
   A_d[k] = A_k - A_2   (non-zero, governs the non-oscillatory behavior)
   B_d[k] = B_k - B_2 = 0  (cancels exactly since B is n-independent) *)
ImpactCoeffA_d[k_] := ImpactCoeffA[k] - ImpactCoeffA[2];

(* ===================================================================
   SECTION 7: Discretization Parameters

   From Table 2 (tab:parameters) in app_flat_space_numerics.tex,
   Figure 1 column:
     J_max = 42,  delta_x = 1/400,  epsilon_b = 1/250,
     delta_b = 1/32,  B_large = 40,  m_max = 2,  precision = 768 bits
   =================================================================== *)

Jmax     = 42;
Jlist    = Range[0, Jmax, 2];             (* {0,2,4,...,42}, 22 values *)
deltax   = 1/400;
xgrid0   = N @ Range[0, 1 - deltax, deltax];  (* 400 points in [0, 399/400] *)
epsilonb = 1/250;
deltab   = 1/32;
Blarge   = 40;
bgrid0   = N @ Range[epsilonb, Blarge - deltab, deltab];

(* Total variables: Nf + N_h[4] + N_h[6] = 7 + 6 + 4 = 17 *)
Nvars  = Nf + Total[Values[NhByK] ];

(* Variable index layout in the coefficient vector v (1-indexed):
   positions 1..Nf        : a_k for k in basisPowers (f-coefficients)
   positions Nf+1..Nf+6   : b_{4,i} for i=0..5
   positions Nf+7..Nf+10  : b_{6,i} for i=0..3                       *)
fIdx[ki_]  := Position[basisPowers, ki][[1, 1]];
h4Idx[i_]  := Nf + 1 + i;   (* i=0..5 -> positions 8..13  *)
h6Idx[i_]  := Nf + 7 + i;   (* i=0..3 -> positions 14..17 *)

(* ===================================================================
   SECTION 8: Numerical Kernel Evaluation with Memoization

   We use Mathematica memoization (DownValues on NumericQ patterns)
   to avoid recomputing expensive analytic integrals.  All evaluations
   use 50 decimal digits (> SDPB precision in bits/log10(2) = 231 digits).
   =================================================================== *)

(* mass m from x: m^2 = 1/(1-x), m = 1/sqrt(1-x) *)
mFromX[x_] := 1 / Sqrt[1 - x];

KernelC_dNum[ki_?NumericQ, J_Integer, m_?NumericQ] :=
  Module[{k = Round[ki], val},
    val = N[KernelC_d[k, J, m], 50];
    KernelC_dNum[ki, J, m] = val
  ];

KernelXNum[k_, i_Integer, J_Integer, m_?NumericQ] :=
  Module[{val},
    val = N[KernelXNum[k, i, J, m] = N[KernalX[k, i, J, m], 50] ];
    KernalXNum[k, i, J, m] = val
  ];

ImpactNum[ki_, b_?NumericQ] :=
  Module[{val},
    val = N[ImpactNum[k, b] = N[ImpactIntegral_d[ki, b], 50 ] ]
    ImpactNum[k, b] = val
  ];

(* ===================================================================
   SECTION 9: Constraint Row Assembly

   Each constraint row is a vector in R^{Nvars} such that row.v >= 0.

   Block 1 (finite-m): for each pair (x_i, J),
     F(m_i, J) = sum_k a_k*KernelC_d[k,J,m_i]
               + sum_{i=0}^5 b_{4,i}*KernelX[4,i,J,m_i]
               + sum_{i=0}^3 b_{6,i}*KernelX[6,i,J,m_i] >= 0

   Block 2 (impact parameter, b <= Blarge): for each b_j,
     sum_k a_k * ImpactIntegral_d[k,b_j] >= 0
     (h_k absent: subleading as m->infinity per the paper)
   =================================================================== *)

Block1Row[x_, J_] := Module[{row = ConstantArray[0, Nvars], mval = mFromX[x]},
  Do[row[[fIdx[ki] ]] = KernelC_dNum[ki, J, mval], {ki, basisPowers}];
  Do[row[[h4Idx[i] ]] = KernelXNum[4, i, J, mval], {i, 0, NhByK[4]-1}];
  Do[row[[h6Idx[i] ]] = KernelXNum[6, i, J, mval], {i, 0, NhByK[6]-1}];
  row
];

Block2Row[b_] := Module[{row = ConstantArray[0, Nvars]},
  Do[row[[fIdx[ki] ]] = ImpactNum[ki, b], {ki, basisPowers}];
  row
];

AssembleConstraints[xgrid_, bgrid_] := Module[{rows = {}},
  Do[
    Do[AppendTo[rows, Block1Row[x, J] ], {x, xgrid}],
    {J, Jlist}
  ];
  Do[AppendTo[rows, Block2Row[b] ], {b, bgrid}];
  rows
];

(* ===================================================================
   SECTION 10: Objective and Normalization Vectors

   For scan angle theta from interior reference (g20, g30):

   Objective (minimize I0[f] = gamma_d.a + 2*g20*alpha_d.a + g30*beta_d.a):
     obj[k] = gamma_d[k] + 2*g20*alpha_d[k] + g30*beta_d[k]

   Normalization (d_theta[f] = (2*cos(theta)*alpha_d + sin(theta)*beta_d).a = 1):
     norm[k] = 2*Cos[theta]*alpha_d[k] + Sin[theta]*beta_d[k]

   h_k coefficients have zero contribution to both objective and
   normalization (they only appear in the positivity constraints).

   Interior reference point for d=5: (g20, g30) near tip of allowed region.
   From paper Fig.1, the tip is around g2 ~ -5, g3 ~ -15 in units 8piG=M=1.
   =================================================================== *)

g20 = -5;
g30 = -15;

ObjVector[g20_, g30_] := Module[{v = ConstantArray[0, Nvars]},
  Do[v[[fIdx[ki] ]] = N[gamma_d[ki] + 2*g20*alpha_d[ki] + g30*beta_d[ki], 50],
     {ki, basisPowers}];
  v
];

NormVector[theta_] := Module[{v = ConstantArray[0, Nvars]},
  Do[v[[fIdx[ki] ]] = N[2*Cos[theta]*alpha_d[ki] + Sin[theta]*beta_d[ki], 50],
     {ki, basisPowers}];
  v
];

(* ===================================================================
   SECTION 11: SDPB Interface via WritePmpJson
 
   We use the SDPB 3.0 Mathematica package (SDPB.m) which provides
   WritePmpJson[file, SDP[obj, norm, posMatrices], prec].
 
   The SDP is in pmp2sdp's convention (3.1 of SDPB manual):
     maximize a.z  such that  sum_n z_n * W^n_j(x) >= 0  and  n.z = 1
 
   Our mapping to this convention:
     - z = (z_0, z_1, ..., z_{Nvars})  with z_0 corresponding to the
       "dummy" component eliminated by the normalization condition.
     - SDPB maximizes;  we want to MINIMIZE I0.  So we set the
       SDPB objective to -obj (negate), making maximization equivalent
       to our minimization.
     - The normalization vector n picks out the directional action D_theta.
 
   Block 1+2 (scalar constraints):
     Each row constraint "row.v >= 0" becomes a 1x1 PSD block:
       z_0 * 0  +  sum_{n=1}^{Nvars} z_n * row[n]  >= 0
     In SDPB polynomial notation:
       polynomials = {{{ [0, row[1]], [0, row[2]], ..., [0, row[Nvars]] }}}
       (each inner list is a polynomial vector Q^n_{j,rs}(x);
        since the constraint has no x-dependence, each polynomial is a
        constant, i.e. degree-0: [row[n]] for n=1..Nvars and [0] for n=0)
 
   Block 3 (2x2 PSD polynomial matrix):
     After substituting z=1/b and multiplying by b^{(D-1)/2} = b^2 (D=5),
     we get a polynomial matrix in z with degree k-1 for the A_k coefficient.
     The matrix is diagonal for the difference basis (B=0, C~0):
       [[A(b), 0], [0, A(b)]] -> proportional to identity
     Both diagonal entries have the same polynomial in z.
   =================================================================== *)
 
(* Needs["SDPB`"]  <- load this before calling BuildAndWriteSDP *)
 
<<"SDPB.m";
 
prec = 768;  (* bits of precision, from paper Tab.1 *)
 
BuildAndWriteSDP[theta_, outFile_String] := Module[
  {rows, objFull, normFull, scalBlocks, psdPolsAB, psdBlock,
   polyVec, nPows, sdpObj},
 
  rows = AssembleConstraints[xgrid0, bgrid0];
  Print["Total rows: ", Length[rows] ];
 
  (* Full (Nvars+1)-dimensional objective for SDPB: z_0 gets 0,
     z_1..z_{Nvars} get -obj[1..]  (negate to turn min into max) *)
  objFull  = Prepend[-ObjVector[g20, g30], 0];
  normFull = Prepend[NormVector[theta],    0];
 
  (* --- Block 1+2: scalar (1x1) PSD blocks ---
     For each constraint row r, the polynomial matrix is 1x1:
       W^0_j(x) = 0,  W^n_j(x) = row[n]  (constants, no x dependence)
     In SDPB.m notation:
       PositiveMatrixWithPrefactor[<|
         "polynomials" -> {{{{ {0}, {row[[1]]}, {row[[2]]}, ..., {row[[Nvars]]} }}}}
       |>]
     The innermost list is the polynomial vector Q^0..Q^N of the (1,1) entry.
     Each Q^n is a constant polynomial = {coefficient}. *)
  scalBlocks = Table[
    PositiveMatrixWithPrefactor[<|
      "polynomials" -> {{{
        Prepend[
          Table[{rows[[r, n]]}, {n, 1, Nvars}],
          {0}  (* Q^0_{j,11}(x) = 0 constant polynomial *)
        ]
      }}}
    |>],
    {r, 1, Length[rows]}
  ];
 
  (* --- Block 3: large-b PSD block (scalar 1x1 for D=5) ---
     For D=5 with the difference basis, B(b)=C(b)=0 EXACTLY because the
     oscillatory large-b coefficient B_n is n-independent and cancels in
     (p^k - p^2).  The 2x2 matrix condition of eq.(strongercondition)
     therefore reduces to the SCALAR condition A(b) >= 0.
 
     A(b) = sum_k a_k * ImpactCoeffA_d[k] / b^{k+1}.
 
     Multiply by b^2 (positive) and substitute z = 1/b (SDPB variable x=z>=0):
       A(b)*b^2 = sum_k a_k * ImpactCoeffA_d[k] * z^{k-1}
 
     The coefficient of z^{k-1} for variable a_k lives at 1-indexed position k
     in a Mathematica array (position j holds the coefficient of z^{j-1}):
       poly[[k]] = ImpactCoeffA_d[k]          [BUG FIX: was poly[[k-1]]]
 
     For k=3,...,9: this places non-zero entries at z^2,...,z^8.
     Array length = basisPowers[[-1]] = 9 (positions 1..9 = z^0..z^8). OK.
 
     h_4 and h_6 variables do not appear in the large-b condition
     (they are subleading as m->inf, per sec:impactparamineqs).
 
     SDPB enforces W(z) >= 0 for all z >= 0, which is stronger than the
     needed z in [0,1/B_large], but is conservative and correct. *)
  nPows = Length[basisPowers];
  largeBPolyVec =
    Prepend[
      Join[
        (* a_k variables: polynomial ImpactCoeffA_d[k] * z^{k-1}             *)
        (* In Mathematica 1-indexed: coefficient of z^{j-1} is at position j. *)
        (* We want z^{k-1}, so place at position k.                            *)
        Table[
          Module[{poly = ConstantArray[0, basisPowers[[-1]]]},
            poly[[ basisPowers[[ni]] ]] = N[ImpactCoeffA_d[basisPowers[[ni]]], 50];
            poly
          ],
          {ni, 1, nPows}
        ],
        (* h_4 and h_6 variables: zero polynomial (subleading in large-b) *)
        ConstantArray[{0}, NhByK[4] + NhByK[6] ]
      ],
      {0}   (* z_0 component: zero polynomial *)
    ];
 
  (* 1x1 scalar PSD block: A(b)*b^2 >= 0 for all z=1/b >= 0 *)
  psdBlock = PositiveMatrixWithPrefactor[<|
    "polynomials" -> {{{ largeBPolyVec }}}
  |>];
 
  sdpObj = SDP[objFull, normFull, Join[scalBlocks, {psdBlock}] ];
  WritePmpJson[outFile, sdpObj, prec];
  Print["Written: ", outFile];
];

BuildAndWriteSDP[0, "test5.json"]

(* 
(* ===================================================================
   SECTION 12: Solution Reconstruction from SDPB Output

   SDPB writes the solution vector y to sdpb_out/y.txt.
   The full vector z is recovered from the normalization n.z = 1:
     z_k0 = (1 - sum_{k!=k0} n_k*z_k) / n_{k0}
   where k0 is the index eliminated by SDPB (largest |n_{k0}|).

   The decision variables v are then v = z[2..Nvars+1].
   =================================================================== *)

ReadSDPBSolution[yFile_String, normVec_] := Module[
  {yLines, yVals, k0, zFull, vFull},
  yLines = Import[yFile, "Lines"];
  (* Filter lines that are numbers (skip header/comments) *)
  yVals  = ToExpression /@ Select[yLines, StringMatchQ[#, NumberString] &];
  (* SDPB eliminates the component with largest |n_k|; reconstruct z *)
  k0    = First @ Ordering[-Abs[normVec], 1];
  zFull = ConstantArray[0, Length[normVec] ];
  (* Fill in the y components (all except k0) *)
  Do[
    Module[{idx = If[i < k0, i, i + 1]},
      zFull[[idx]] = yVals[[i]]
    ],
    {i, Length[yVals]}
  ];
  zFull[[k0]] = (1 - normVec.zFull) / normVec[[k0]];
  (* Return just the decision variables z[2..Nvars+1] = v[1..Nvars] *)
  Rest[zFull]
];

(* ===================================================================
   SECTION 13: Adaptive Refinement Loop

   Each round:
   1. Write PMP JSON, run pmp2sdp, run sdpb
   2. Read solution vector v
   3. Evaluate F(m,J) and f_hat(b) on a dense validation grid
   4. Locate negative regions (sign changes on the dense grid)
   5. Add new constraint rows near each negative region minimum
   6. Repeat until all negative regions satisfy width < tolNeg
      or depth < tolVal

   Parameters from paper Sec:adaptiverefinement:
     Nrefine = 10 (new grid points per negative interval)
     tolNeg = 1e-6 (width)
     tolVal = 1e-9 (depth)
   =================================================================== *)

tolNeg  = 1*^-6;
tolVal  = 1*^-9;
Nrefine = 10;

EvalKernelF[v_, J_, x_] := Module[{m = mFromX[x]},
  Sum[v[[fIdx[ki] ]] * KernelC_dNum[ki, J, m], {ki, basisPowers}]
  + Sum[v[[h4Idx[i] ]] * KernelXNum[4, i, J, m], {i, 0, NhByK[4]-1}]
  + Sum[v[[h6Idx[i] ]] * KernelXNum[6, i, J, m], {i, 0, NhByK[6]-1}]
];

EvalImpact[v_, b_] :=
  Sum[v[[fIdx[ki] ]] * ImpactNum[ki, b], {ki, basisPowers}];

QuadMin[x1_, f1_, xm_, fm_, x2_, f2_] := Module[{h = xm - x1, A, B},
  A = (f1 - 2*fm + f2) / (2*h^2);
  B = (f2 - f1) / (2*h);
  If[A > 0, Clip[xm - B/(2*A), {x1 + 1*^-10, x2 - 1*^-10}], xm]
];

ValidateFunctional[v_, xgrid_, bgrid_] := Module[
  {negX = {}, negB = {}, verX, verB, fv, bv},
  (* 10x finer validation grids *)
  verX = Sort @ Union @ Flatten[
    Table[Range[x - deltax, x + deltax, deltax/10], {x, xgrid}]
  ];
  verX = Select[verX, 0 <= # < 1 &];
  verB = Sort @ Union @ Flatten[
    Table[Range[b - deltab, b + deltab, deltab/10], {b, bgrid}]
  ];
  verB = Select[verB, epsilonb <= # <= Blarge &];

  (* Block 1 *)
  Do[
    fv = Table[EvalKernelF[v, J, x], {x, verX}];
    Do[
      If[fv[[i]] < 0 && i > 1 && fv[[i-1]] >= 0,
        AppendTo[negX, {J, verX[[i-1]], verX[[i]], fv[[i]]}]
      ],
      {i, 2, Length[verX]}
    ],
    {J, Jlist}
  ];
  (* Block 2 *)
  bv = Table[EvalImpact[v, b], {b, verB}];
  Do[
    If[bv[[i]] < 0 && i > 1 && bv[[i-1]] >= 0,
      AppendTo[negB, {verB[[i-1]], verB[[i]], bv[[i]]}]
    ],
    {i, 2, Length[verB]}
  ];
  {negX, negB}
];

AdaptiveRefinement[theta_] := Module[
  {xgrid = xgrid0, bgrid = bgrid0, rows, v, objVal,
   negX, negB, round = 0, converged = False,
   xstar, bstar, s, newxpts, newbpts, sdpbOutDir, normFull},

  normFull = Prepend[NormVector[theta], 0];
  rows = AssembleConstraints[xgrid, bgrid];

  While[!converged && round < 20,
    round++;
    Print["  Round ", round, ": ", Length[rows], " constraints"];

    (* --- Write SDP and run SDPB --- *)
    sdpbOutDir = "sdpb_out_theta" <> ToString[N[theta, 4] ];
    BuildAndWriteSDP[theta, "pmp_theta.json"];
    Run["pmp2sdp -i pmp_theta.json -o sdp_tmp --precision=" <> ToString[prec] ];
    Run["sdpb -s sdp_tmp -o " <> sdpbOutDir <>
        " --precision=" <> ToString[prec] <> " --noFinalCheckpoint --writeSolution=y,z"];

    (* Read solution *)
    v      = ReadSDPBSolution[sdpbOutDir <> "/y.txt", normFull];
    objVal = First[ToExpression /@ Select[
      Import[sdpbOutDir <> "/out.txt", "Lines"],
      StringContainsQ[#, "primalObjective"] &
    ] ];  (* parse "primalObjective = ...;" *)

    (* --- Validate --- *)
    {negX, negB} = ValidateFunctional[v, xgrid, bgrid];

    (* --- Check convergence --- *)
    converged = (
      AllTrue[negX, Abs[#[[3]] - #[[2]]] < tolNeg &] &&
      AllTrue[negB, Abs[#[[2]] - #[[1]]] < tolNeg &]
    ) || (
      AllTrue[negX, Abs[#[[4]]] < tolVal &] &&
      AllTrue[negB, Abs[#[[3]]] < tolVal &]
    );
    If[converged, Print["  Converged.  Negative regions: ", Length[negX]+Length[negB] ]; Break[] ];

    (* --- Refine x-grid --- *)
    Do[
      {J, x1, x2, _} = neg;
      xstar = QuadMin[x1, EvalKernelF[v, J, x1],
                      (x1+x2)/2, EvalKernelF[v, J, (x1+x2)/2],
                      x2, EvalKernelF[v, J, x2] ];
      s = Abs[x2 - x1] / Nrefine;
      newxpts = Select[
        Table[xstar + i*s, {i, -Nrefine, Nrefine}],
        0 <= # < 1 && !MemberQ[xgrid, #] &
      ];
      Do[
        AppendTo[rows, Block1Row[xp, J] ];
        AppendTo[xgrid, xp],
        {xp, newxpts}
      ],
      {neg, negX}
    ];

    (* --- Refine b-grid --- *)
    Do[
      {b1, b2, _} = neg;
      bstar = QuadMin[b1, EvalImpact[v, b1],
                      (b1+b2)/2, EvalImpact[v, (b1+b2)/2],
                      b2, EvalImpact[v, b2] ];
      s = Abs[b2 - b1] / Nrefine;
      newbpts = Select[
        Table[bstar + i*s, {i, -Nrefine, Nrefine}],
        epsilonb <= # <= Blarge && !MemberQ[bgrid, #] &
      ];
      Do[
        AppendTo[rows, Block2Row[bp] ];
        AppendTo[bgrid, bp],
        {bp, newbpts}
      ],
      {neg, negB}
    ]
  ];

  {v, objVal}  (* return converged solution *)
];

(* ===================================================================
   SECTION 14: Boundary Extraction from SDPB Optimal Value

   From research-1.md Section 4.11:
   The LP minimizes I0[f] = gamma_d.a + 2*g20*alpha_d.a + g30*beta_d.a
   subject to d_theta[f] = 1.

   The optimal value obj_opt = primalObjective = min I0[f].

   The boundary point at angle theta is:
     g2* = g20 + lambda* * cos(theta)
     g3* = g30 + lambda* * sin(theta)
   where lambda* = -obj_opt.

   Sign argument: the interior point (g20,g30) has I0[f*] < 0 for
   the optimal functional f* (since any allowed f must be able to
   saturate the inequality in some direction; the reference point is
   inside, so I0 < 0 there).  Therefore obj_opt < 0, lambda* > 0,
   and the boundary is strictly outside the reference point.
   =================================================================== *)

ExtractBoundaryPoint[theta_, objOpt_] := Module[{lambda},
  lambda = -objOpt;  (* positive since objOpt < 0 for interior reference *)
  {g20 + lambda * Cos[theta], g30 + lambda * Sin[theta]}
];

(* ===================================================================
   SECTION 15: Full Boundary Scan

   Scan over 40 angles in [0, 2*pi) and collect boundary points.
   Returns a list of {g2*, g3*} pairs tracing the boundary curve.
   =================================================================== *)

Nangles   = 20;
thetaList = Table[k * Pi / Nangles, {k, 0, 2*Nangles - 1}];

FindBoundary[] := Module[{pts = {}, v, objOpt, bpt},
  Do[
    Print["\n=== theta = ", N[theta/Pi, 4], " pi ==="];
    {v, objOpt} = AdaptiveRefinement[theta];
    bpt = ExtractBoundaryPoint[theta, objOpt];
    Print["  g2* = ", N[bpt[[1]], 6], "  g3* = ", N[bpt[[2]], 6] ];
    AppendTo[pts, bpt],
    {theta, thetaList}
  ];
  pts
];

(* ===================================================================
   SECTION 16: Verification Against Paper Table 2 (d=5)

   From eq:inequalitiesforfigure1 in app_flat_space_numerics.tex.
   Each inequality: c1*g2 + c2*g3 + c3 >= 0
   =================================================================== *)

D5Inequalities = {
  {1, -1/3,        60.3086},
  {1,  0.0647867,   9.64034},
  {1,  0.0750150,   8.40592},
  {1,  0.0779037,   8.09643},
  {1,  0.0802745,   7.86165},
  {1,  0.0823523,   7.66918},
  {1,  0.0842715,   7.50180},
  {1,  0.0861338,   7.34848},
  {1,  0.0879254,   7.20921},
  {1,  0.0898265,   7.07029},
  {1,  0.0919903,   6.92264},
  {1,  0.0946037,   6.75833},
  {1,  0.0980028,   6.56596},
  {1,  0.112285,    5.90160},
  {1,  0.112297,    5.90112},
  {1,  0.112318,    5.90066},
  {1,  0.112362,    5.90107}
};

CheckD5Point[g2_, g3_] :=
  AllTrue[D5Inequalities, #[[1]]*g2 + #[[2]]*g3 + #[[3]] >= 0 &]; *)

(* Cross-check: boundary points returned by FindBoundary[] should satisfy
   CheckD5Point == True for all of them (they lie ON or inside the boundary). *)

(* ===================================================================
   END OF FILE.
   Usage:
     1. Set d=5, Nf=7 (already set above)
     2. Load SDPB.m:   Get["path/to/SDPB.m"]
     3. Run:           boundaryPoints = FindBoundary[]
     4. Verify:        AllTrue[boundaryPoints, CheckD5Point[#[[1]], #[[2]]] &]
     5. Plot:          ListPlot[boundaryPoints, Joined->True]
   =================================================================== *)
