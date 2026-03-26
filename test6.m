(* ===========================================================
   SYMBOLIC PHASE DESIGN PRINCIPLES
   -----------------------------------------------------------
   Problem: FullSimplify recognizes polynomial patterns and
   rewrites them as GegenbauerC[n, λ, expr], reintroducing
   the very head we are trying to eliminate. Once GegenbauerC
   survives into xkExpr, N[..., 60] converts its exact integer
   parameters (2, 1) to 60-digit floats (2.000..., 1.000...),
   and Mathematica has no numerical evaluation rule for
   GegenbauerC^(0,0,k)[2.000..., 1.000..., 1.000...], so
   NIntegrate sees non-numerical values.

   Solution: use Together (rational arithmetic only) instead
   of FullSimplify for the cached expression. Together never
   recognizes special-function patterns.

   PJ must return an explicit polynomial in its argument
   with exact rational coefficients — never a GegenbauerC
   head — so that Residue differentiates plain rational
   functions, not special functions.
   =========================================================== *)

ClearAll["Global`*"]
Remove["Global`*"]

d  = 5;
Nf = 7;

(* PJ[J, x]: explicit polynomial ratio with exact integer
   (d-3)/2 = 1, so GegenbauerC[J, 1, x] is a polynomial in x.
   Mathematica evaluates GegenbauerC[n, integer, x] to an
   explicit polynomial for concrete integer n.
   We call Expand to be sure the polynomial form is manifest,
   and Simplify (not FullSimplify) on the ratio. *)
PJ[J_Integer, x_] :=
  Expand[GegenbauerC[J, (d-3)/2, x]] /
  Expand[GegenbauerC[J, (d-3)/2, 1]];

(* Sanity check: PJ should return a pure polynomial in x,
   never a GegenbauerC head. *)
With[{testExpr = PJ[2, testX]},
  If[! FreeQ[testExpr, GegenbauerC],
    Print["ERROR: PJ still contains GegenbauerC — fix PJ definition"],
    Print["OK: PJ[2, x] = ", testExpr]
  ]
];

PJprime1[J_] := J*(J + d - 3) / (d - 2);


(* ===========================================================
   NUMERIC KERNELS T1, T2, T3
   PJ is called with a purely numeric argument here
   (1 - 2p^2/mP^2), so no symbolic ambiguity arises.
   =========================================================== *)
T1[n_Integer, J_Integer, m_?NumericQ] /; m > 1 && n > 1 :=
  Module[{mP = SetPrecision[m, 60]},
    NIntegrate[
      p^n * (2*mP^2 - p^2) * PJ[J, 1 - 2*p^2/mP^2] /
      (mP^2 * (mP^2 - p^2)^2),
      {p, 0, 1},
      WorkingPrecision -> 60, PrecisionGoal -> 50, MaxRecursion -> 20
    ]
  ];

T2[n_Integer, m_?NumericQ] /; m > 1 && n > 1 :=
  Module[{mP = SetPrecision[m, 60]},
    NIntegrate[
      -p^n * p^4 * (4*mP^2 - 3*p^2) /
      (mP^6 * (mP^2 - p^2)^2),
      {p, 0, 1},
      WorkingPrecision -> 60, PrecisionGoal -> 50, MaxRecursion -> 20
    ]
  ];

T3[n_Integer, J_Integer, m_?NumericQ] /; m > 1 && n > 1 :=
  Module[{mP = SetPrecision[m, 60]},
    PJprime1[J] * NIntegrate[
      4 * p^n * p^6 /
      (mP^6 * (mP^4 - p^4)),
      {p, 0, 1},
      WorkingPrecision -> 60, PrecisionGoal -> 50, MaxRecursion -> 20
    ]
  ];

KernelC[n_, J_, m_] := T1[n, J, m] + T2[n, m] + T3[n, J, m];


gammaF[n_] := 1/(n - 1);
alphaF[n_] := 1/(n + 1);
betaF[n_]  := 1/(n + 3);

basisPowers = Range[3, Nf + 2];    (* {3,4,...,9} for Nf=7 *)

gammaD[k_] := gammaF[k] - gammaF[2];
alphaD[k_] := alphaF[k] - alphaF[2];
betaD[k_]  := betaF[k]  - betaF[2];

KernelCD[k_, J_, m_] := KernelC[k, J, m] - KernelC[2, J, m];


NhByK = <|4 -> 6, 6 -> 4|>;


(* ===========================================================
   XKERNEL COMPONENTS
   Since PJ returns explicit polynomials, these are rational
   functions in their arguments from the start.
   =========================================================== *)
XKernelFirstTerm[k_Integer, m_, J_, u_] :=
  (2*m^2 + u) * m^2 * PJ[J, 1 + 2*u/m^2] /
  (u * m^2 * (m^2 + u))^(1 + k/2);

(* Corrected integrand for the residue term *)
XKernelResidueIntegrand[k_Integer, m_, J_, u_, up_] :=
  ((2*m^2 + up)*(m^2 - up)*(m^2 + 2*up)) /
  ((u - up)*up*(m^2 - u)*(m^2 + up)*(m^2 + u + up)) *
  PJ[J, 1 + 2*up/m^2] / (up*m^2*(m^2 + up))^(k/2);

XKernelCached[k_Integer, J_Integer] :=
  XKernelCached[k, J] =
    Module[{ft, ri, ser, res, expr},
      ft = XKernelFirstTerm[k, mSym, J, uSym];
      ri = XKernelResidueIntegrand[k, mSym, J, uSym, up];
      (* Compute residue using series expansion *)
      ser = Series[ri, {up, 0, -1}];
      res = Coefficient[ser, up^(-1)];
      expr = Together[ft - res];
      (* Optional: check analyticity at u=0 *)
      If[! FreeQ[Series[expr, {uSym, 0, 0}], uSym^(-1)],
        Print["Warning: expression still has a pole at u=0."];
      ];
      If[! FreeQ[expr, GegenbauerC],
        Print["FATAL: GegenbauerC survived."]; Abort[]
      ];
      expr
    ];

(* --- DEBUG: print the expression and its series at u=0 --- *)
Print["XKernelCached[4,4] = ", XKernelCached[4,4]];
Print["Series around u=0: ", Series[XKernelCached[4,4], {uSym, 0, 2}]];

(* Run the cache-build test and verify both conditions *)
With[{testExpr = XKernelCached[4, 4]},
  Print["FreeQ[Residue]:     ", FreeQ[testExpr, Residue]];
  Print["FreeQ[GegenbauerC]: ", FreeQ[testExpr, GegenbauerC]];
];


(* ===========================================================
   KERNELXNUM: NUMERIC EVALUATION
   -----------------------------------------------------------
   xkExpr is a pure rational function of mSym and uSym
   (no GegenbauerC anywhere, verified above).

   Substitution order:
     1. Substitute mSym -> mP, uSym -> -pp^2 symbolically.
        This keeps exact rational arithmetic as long as mP is
        a SetPrecision number and pp is still a formal symbol.
     2. Wrap in Function[pp, Evaluate[...]] so the rational
        expression is compiled once with mP baked in, and pp
        is the only remaining variable. Evaluate[] fires at
        Function-construction time (mP is already bound),
        not at NIntegrate sampling time.
     3. NIntegrate calls integrandFn[numericPP]: substitutes
        a 60-digit float for pp into a closed rational
        expression — fully numerical, no special functions.

   We do NOT call N[..., 60] around the whole substitution
   because that would convert exact integer literal arguments
   of any surviving GegenbauerC (if our assertion somehow
   failed) into floats, hiding the bug and breaking evaluation.
   The 60-digit precision comes from mP = SetPrecision[m, 60].
   =========================================================== *)
KernelXNum[k_Integer, n_Integer, J_Integer, m_?NumericQ] /; m > 1 :=
  KernelXNum[k, n, J, m] =
  Module[{mP = SetPrecision[m, 60], xkExpr, integrandFn, eps = 10^-10},
    xkExpr = XKernelCached[k, J];
    integrandFn = Function[pp,
      Evaluate[pp^n * (xkExpr /. {mSym -> mP, uSym -> -pp^2})]
    ];
    NIntegrate[
      integrandFn[pp],
      {pp, eps, 1},   (* avoid pp=0 *)
      WorkingPrecision -> 60, PrecisionGoal -> 50, MaxRecursion -> 30,
      Method -> {"GlobalAdaptive", "SingularityHandler" -> None}
    ]
  ];


(* ===========================================================
   IMPACT INTEGRALS (unchanged)
   =========================================================== *)
ImpactIntegral[n_, b_] :=
  HypergeometricPFQ[{(n+1)/2}, {(d-2)/2, (n+3)/2}, -b^2/4] / (n+1);

ImpactIntegralD[n_, b_] :=
  ImpactIntegral[n, b] - ImpactIntegral[2, b];

ImpactCoeffA[n_] :=
  2^n * Gamma[(d-2)/2] * Gamma[(n+1)/2] / Gamma[(d-n-3)/2];

ImpactCoeffB[n_] :=
  2^((d-3)/2) * Gamma[(d-2)/2] / Sqrt[Pi];

ImpactCoeffAD[k_] := ImpactCoeffA[k] - ImpactCoeffA[2];


(* ===========================================================
   GRID AND INDEX SETUP (unchanged)
   =========================================================== *)
(* Jmax     = 42; *)
Jmax = 0;
Jlist    = Range[0, Jmax, 2];
(* deltax   = 1/400; *)
deltax = 1/10;
xgrid0   = N @ Range[0, 1 - deltax, deltax];
epsilonb = 1/250;
deltab   = 1/32;
Blarge   = 40;
bgrid0   = N @ Range[epsilonb, Blarge - deltab, deltab];

Nvars = Nf + Total[Values[NhByK]];

fIdx[ki_] := Position[basisPowers, ki][[1, 1]];
h4Idx[i_] := Nf + 1 + i;
h6Idx[i_] := Nf + 7 + i;

mFromX[x_] := 1 / Sqrt[1 - x];


(* ===========================================================
   CACHED NUMERIC KERNELS (unchanged logic, precision restored)
   =========================================================== *)
KernelCdNum[ki_, J_, m_?NumericQ] :=
  Module[{val},
    val = N[KernelCD[ki, J, m], 60];
    KernelCdNum[ki, J, m] = val
  ];

ImpactNum[ki_, b_?NumericQ] :=
  Module[{val},
    val = N[ImpactIntegralD[ki, b], 60];
    ImpactNum[ki, b] = val
  ];


(* ===========================================================
   ROW ASSEMBLERS AND CONSTRAINT MATRIX (unchanged)
   =========================================================== *)
Block1Row[x_, J_] :=
  Module[{row = ConstantArray[0, Nvars], mval = mFromX[x]},
    Do[row[[fIdx[ki]]] = KernelCdNum[ki, J, mval],  {ki, basisPowers}];
    Do[row[[h4Idx[i]]] = KernelXNum[4, i, J, mval], {i, 0, NhByK[4]-1}];
    Do[row[[h6Idx[i]]] = KernelXNum[6, i, J, mval], {i, 0, NhByK[6]-1}];
    row
  ];

Block2Row[b_] :=
  Module[{row = ConstantArray[0, Nvars]},
    Do[row[[fIdx[ki]]] = ImpactNum[ki, b], {ki, basisPowers}];
    row
  ];

AssembleConstraints[xgrid_, bgrid_] :=
  Module[{rows = {}},
    Do[
      Do[AppendTo[rows, Block1Row[x, J]], {x, xgrid}],
      {J, Jlist}
    ];
    Do[AppendTo[rows, Block2Row[b]], {b, bgrid}];
    rows
  ];


(* ===========================================================
   OBJECTIVE AND NORM VECTORS (unchanged)
   =========================================================== *)
g20 = -5;
g30 = -15;

ObjVector[g20_, g30_] :=
  Module[{v = ConstantArray[0, Nvars]},
    Do[v[[fIdx[ki]]] = N[gammaD[ki] + 2*g20*alphaD[ki] + g30*betaD[ki], 60],
       {ki, basisPowers}];
    v
  ];

NormVector[theta_] :=
  Module[{v = ConstantArray[0, Nvars]},
    Do[v[[fIdx[ki]]] = N[2*Cos[theta]*alphaD[ki] + Sin[theta]*betaD[ki], 60],
       {ki, basisPowers}];
    v
  ];


(* ===========================================================
   SMOKE TEST — runs before the heavy SDPB machinery
   =========================================================== *)
(* Print["Testing KernelXNum[4, 1, 2, 1.5] ..."];
result = KernelXNum[4, 1, 4, 1.1];
Print["result = ", result];

Print["Testing KernelCdNum[ki_, J_, m_?NumericQ]..."];
result = KernelCdNum[4, 2, 1.1];
Print["result = ", result];

Print["Testing ImpactNum[ki_, b_?NumericQ]..."];
result = ImpactNum[4, 0.5];
Print["result = ", result]; *)

(* ===========================================================
   SDPB OUTPUT (unchanged)
   =========================================================== *)
<<"SDPB.m";

prec = 768;

BuildAndWriteSDP[theta_, outFile_String, xgrid_, bgrid_] :=
  Module[{rows, objFull, normFull, scalBlocks, psdBlock,
          polyVec, nPows, sdpObj},

    rows = AssembleConstraints[xgrid0, bgrid0];
    Print["Total rows: ", Length[rows]];

    objFull  = Prepend[-ObjVector[g20, g30], 0];
    normFull = Prepend[-NormVector[theta],   0];

    scalBlocks = Table[
      PositiveMatrixWithPrefactor[<|
        "polynomials" -> {{{
          Prepend[
            Table[{rows[[r, n]]}, {n, 1, Nvars}],
            {0}
          ]
        }}}
      |>],
      {r, 1, Length[rows]}
    ];

    nPows   = Length[basisPowers];
    polyVec = Function[{matEntry},
      Prepend[
        Join[
          Table[
            Module[{poly = ConstantArray[0, basisPowers[[-1]]]},
              poly[[basisPowers[[ni]]]] = N[ImpactCoeffAD[basisPowers[[ni]]], 60];
              poly
            ],
            {ni, 1, nPows}
          ],
          ConstantArray[{0}, NhByK[4] + NhByK[6]]
        ],
        {0}
      ]
    ];

    psdBlock = PositiveMatrixWithPrefactor[<|
      "polynomials" -> {
        {
          {polyVec["M11"], Prepend[ConstantArray[{0}, Nvars], {0}]},
          {Prepend[ConstantArray[{0}, Nvars], {0}], polyVec["M22"]}
        }
      }
    |>];

    sdpObj = SDP[objFull, normFull, Join[scalBlocks, {psdBlock}]];
    WritePmpJson[outFile, sdpObj, prec];
    Print["Written: ", outFile];
  ];

BuildAndWriteSDP[0.0, "test5.json", xgrid0, bgrid0];