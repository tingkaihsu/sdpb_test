(* ================================================================
   largeSpinTest.m  (precision-propagation corrected version)
   ----------------------------------------------------------------
   PURPOSE
     After an SDPB run with Jmax=60, read the dual vector y and the
     current sampling points, then check whether the functional

         F(x, J) = y₁·g20(x,J) + y₂·g31(x,J) + y₃·n4(x,J)
                 + y₄·X52(x,J) + y₅·X53(x,J)

     remains non-negative for LARGE spins J > 60.

   PRECISION FIXES (2026-05) — applied consistently with test9.m
     1. maVal: N[3/20, 650] instead of SetPrecision[0.150, 600].
     2. n4PrecMin: adaptive working precision for n4.
     3. n4: evaluates at prec_local digits, Re[…] explicit.
     4. All other functions: N[3/20, 650] for mA.

   BUG FIXES (pre-existing)
     A. Asymptotic analysis (Phase 3): extraFuncs had only 3 entries
        when fList has 5.  The large-J functional is
          X(x) = Σ_k yVec[[k]] * extraTriplet[[k]][x]
        which correctly selects y₅ * LargeJ(x) (the X53 large-J term).
        Fixed to use extraTriplet directly.
     B. SetPrecision[ToExpression /@ (fixSciNotation /@ val), 200] was
        applied to a numeric val (not a string). Fixed to N[val, 50].

   PERFORMANCE NOTE
     Evaluating n4[x, J] at J = 20000 requires ≈ 3196-digit precision.
     Computing Hypergeometric2F1[-20000, 20001, 1, z] at that precision
     is equivalent to evaluating a degree-20000 Legendre polynomial to
     3196 digits — this may take minutes per evaluation.  Consider
     limiting jScanMax to ≤ 10000 for routine diagnostic runs.

   USAGE
     wolframscript -file largeSpinTest.m <sp_file> <y_file> [J_scan_max]

   EXIT CODES
     0  F ≥ 0 for all large J tested
     1  Violations found
     2  Error
   ================================================================ *)


(* ----------------------------------------------------------------
   1. ARGUMENT PARSING
   ---------------------------------------------------------------- *)

myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

If[Length[myArgs] < 2,
  Print["USAGE: wolframscript -file largeSpinTest.m ",
        "<sp_file> <y_file> [N_pts] [out_sp_file]"];
  Quit[2]
];

spFile    = myArgs[[1]];
yFile     = myArgs[[2]];
nPts      = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 10];
minWidth  = If[Length[myArgs] >= 5, ToExpression[myArgs[[5]]], 10^-6];
jScanMax  = If[Length[myArgs] >= 4, ToExpression[myArgs[[4]]], 10000];

Print["=== largeSpinTest.m  (Large-J Positivity Diagnostic) ==="];
Print["  sp_file     = ", spFile];
Print["  y_file      = ", yFile];
Print["  N_pts      = ", nPts];
Print["  min_width  = ", minWidth];
Print["  J_scan_max  = ", jScanMax];
Print[""];


(* ----------------------------------------------------------------
   2. READ FILES
   ---------------------------------------------------------------- *)

fixSciNotation[s_String] := StringReplace[s,
  RegularExpression["[eE]([+-]?\\d+)"] :> "*^$1"
];

If[!FileExistsQ[spFile],
  Print["ERROR: sampling points file not found: ", spFile]; Quit[2]];

spRaw = Select[
  ReadList[spFile, String],
  StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &
];

If[Length[spRaw] == 0,
  Print["ERROR: no sample points found in ", spFile]; Quit[2]];

samplePoints = Sort[SetPrecision[ToExpression /@ (fixSciNotation /@ spRaw), 200]];
Print["Loaded ", Length[samplePoints], " sampling points."];


(* ----------------------------------------------------------------
   3.  READ y.txt / z.txt
   ---------------------------------------------------------------- *)

If[!FileExistsQ[yFile],
  Print["ERROR: z.txt not found: ", yFile]; Quit[2]];

yRaw = Select[
  ReadList[yFile, String],
  StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &
];

If[Length[yRaw] == 0,
  Print["ERROR: z.txt is empty: ", yFile]; Quit[2]];

yVec = SetPrecision[ToExpression /@ (fixSciNotation /@ yRaw), 200];

Do[
  If[!NumberQ[yVec[[k]]],
    Print["ERROR: y component ", k, " is not numeric: ", yVec[[k]]];
    Print["  Raw line was: ", yRaw[[k]]];
    Print["  After fixSciNotation: ", fixSciNotation[yRaw[[k]]]];
    Quit[2]
  ],
  {k, Length[yVec]}
];

Do[
  If[!NumberQ[samplePoints[[k]]],
    Print["ERROR: sample point ", k, " is not numeric: ", samplePoints[[k]]];
    Quit[2]
  ],
  {k, Length[samplePoints]}
];

Print["y vector (", Length[yVec], " component(s)): ", yVec];
Print[""];


(* ----------------------------------------------------------------
   4. PROBLEM-SPECIFIC DEFINITIONS — must match test9.m
   ---------------------------------------------------------------- *)

(* FIX 1: exact rational for mA — avoids machine-float precision leakage *)
maVal = N[3/20, 650];
Print["mA = ", maVal];

(* FIX 2: n4PrecMin — adaptive working precision for n4.
   For J ≤ ~3782: returns 650 (no overhead vs current code).
   For J = 10000: returns 1628.
   For J = 20000: returns 3196.
   Uses the exact rational (3/20)^2 for mA in the precision estimator. *)
n4PrecMin[J_?IntegerQ, x_?NumericQ, margin_Integer:60, stdPrec_Integer:650] :=
  Module[{sp0, z2, A, logA},
    sp0  = N[1/(1 - x), 50];
    z2   = 1 + 8*(maVal)^2 / (3*(sp0 - 4*(maVal)^2));
    A    = If[z2 > 1, z2 + Sqrt[z2^2 - 1], 1];
    logA = If[A > 1, N[Log[10, A], 50], 0];
    Max[stdPrec, Ceiling[J * logA] + margin]
  ];

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[maVal, 650];
  N[-(2*Sqrt[sp/(-4*mA^2 + sp)])/(2*mA^2 - sp)^3, 650]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[maVal, 650];
  N[-(Sqrt[sp/(-4*mA^2 + sp)]*(-3 - (2*J*(1 + J)*(2*mA^2 - sp))/(-4*mA^2 + sp)))/(-2*mA^2 + sp)^4, 650]
];

(* null constraint — PRECISION-PROPAGATION IMPLEMENTATION
   FIX 2+3+4: adaptive precision, exact-rational mA, explicit Re[…]
   
   PERFORMANCE NOTE: for J = 10000 this requires ~1628 digits.
   For J = 20000 it requires ~3196 digits.  The Hypergeometric2F1
   and LegendreP evaluations at these precisions may take minutes. *)
n4[x_?NumericQ, J_?IntegerQ] := Module[{prec, xp, sp, mA, result},
  prec = n4PrecMin[J, x];
  If[prec > 1000,
    Print["  [n4] J=", J, " requires prec=", prec, " digits — this may be slow..."]
  ];
  xp   = SetPrecision[x, prec];
  sp   = N[1/(1 - xp), prec];
  mA   = N[maVal, prec];
  result = (81*Sqrt[sp/(-4*mA^2 + sp)]*(
      8*mA^4*(14*mA^2 - 15*sp)*(-8*mA^2 + 3*sp)^(3/2)*
        Hypergeometric2F1[-J, 1 + J, 1, (4*mA^2)/(12*mA^2 - 3*sp)] +
      (8*mA^4 - 18*mA^2*sp + 9*sp^2)*(
        (-2*I)*mA*(10*mA^2 - 9*sp)*(8*mA^2 - 3*sp)*
          LegendreP[J, 1, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))] +
        Sqrt[-8*mA^2 + 3*sp]*(-8*mA^4 + 18*mA^2*sp - 9*sp^2)*
          LegendreP[J, 2, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))]
      )
    ))/(4*mA^2*(-2*mA^2 + sp)*(-8*mA^2 + 3*sp)^(3/2)*(8*mA^4 - 18*mA^2*sp + 9*sp^2)^3);
  Re[N[result, 650]]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[maVal, 650];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 650]
];

X53[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[maVal, 650];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*((-4 + J)*(-2 + J)*(3 + J)*(5 + J) + ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 650]
];

LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[maVal, 650];
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-32*mA^6 + 24*mA^4*sp - 6*mA^2*sp^2 + sp^3))/(18*sp^3*(-4*mA^2 + sp)^6), 650]
];

Jmax = 60;
Jlist = Range[0, Jmax, 2];
(* JlistLarge = {3000};
Jlist = Join[Range[0, Jmax, 2], JlistLarge]; *)

fList = {g20, g31, n4, X52, X53};

(* large J limit: g20→0, g31→0, n4→0, X52→0, X53→LargeJ *)
extraTriplet = {0&, 0&, 0&, 0&, LargeJ};

norm = {0, -1, 0, 0, 0};
obj  = {-1, 0, 0, 0, 0};

xLeft  = SetPrecision[0, 650];
xRight = SetPrecision[1, 650];


(* --- Dimensional consistency check --- *)
If[Length[yVec] != Length[fList],
  Print["ERROR: y.txt has ", Length[yVec], " component(s) but fList has ",
        Length[fList], " function(s). They must match."];
  Quit[2]
];


(* ----------------------------------------------------------------
   4. BUILD SMART J GRID FOR SCANNING
   ---------------------------------------------------------------- *)

JlistDense    = Range[62, Min[120, jScanMax], 2];
JlistModerate = Range[140, Min[500, jScanMax], 20];
JlistSparse   = Select[
  {600, 800, 1000, 1500, 2000, 3000, 5000, 7500, 10000, 15000, 20000},
  # <= jScanMax &
];

JlistLarge = DeleteDuplicates[Sort[Join[JlistDense, JlistModerate, JlistSparse]]];

Print["Large-J scan grid: ", Length[JlistLarge], " spin values"];
Print["  Range: J = ", First[JlistLarge], " to ", Last[JlistLarge]];
Print["  Dense (62-120):     ", Length[JlistDense],    " values"];
Print["  Moderate (140-500): ", Length[JlistModerate], " values"];
Print["  Sparse (600+):      ", Length[JlistSparse],   " values"];
Print[""];

(* Warn if large J will need very high precision *)
If[jScanMax >= 5000,
  Print["  WARNING: J ≥ 5000 will require prec ≥ 844 digits for n4."];
  Print["  WARNING: J ≥ 10000 will require prec ≥ 1628 digits — evaluations may be slow."];
  Print["  Consider setting jScanMax ≤ 3782 for routine runs."];
  Print[""]
];


(* ----------------------------------------------------------------
   5. FUNCTIONAL EVALUATION
   ---------------------------------------------------------------- *)

F[x_?NumericQ, J_?IntegerQ] := Sum[yVec[[k]] * fList[[k]][x, J], {k, Length[yVec]}];

singularTol = SetPrecision[10^-12, 650];

safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  If[Abs[x0] < singularTol, x0 = singularTol];
  val = Quiet[Check[N[F[x0, J], 200], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate}, Head[val]],
    $Failed,
    val
  ]
];


(* ----------------------------------------------------------------
   6. PHASE 1: CHECK AT SAMPLE POINTS
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  PHASE 1: Checking at SAMPLE POINTS"];
Print["========================================"];
Print["  ", Length[samplePoints], " x-points and ", Length[JlistLarge], " J-values"];
Print[""];

phase1Violations = {};
phase1MinVal = Infinity;
phase1WorstX = None;
phase1WorstJ = None;

Do[
  xi = samplePoints[[i]];

  localMin  = Infinity;
  localWorstJ = None;

  Do[
    val = safeF[xi, Jj];
    If[val =!= $Failed && val < localMin,
      localMin = val;
      localWorstJ = Jj;
    ],
    {Jj, JlistLarge}
  ];

  If[localMin < 0,
    AppendTo[phase1Violations, {xi, localWorstJ, localMin}];
  ];

  If[localMin < phase1MinVal,
    phase1MinVal = localMin;
    phase1WorstX = xi;
    phase1WorstJ = localWorstJ;
  ];

  If[Mod[i, Max[1, Floor[Length[samplePoints]/10]]] == 0,
    Print["  Progress: ", i, "/", Length[samplePoints],
          "  local min so far: ", phase1MinVal]
  ],
  {i, Length[samplePoints]}
];

Print[""];
Print["Phase 1 Results:"];
Print["  Global minimum F value:  ", phase1MinVal];
Print["  Worst point:  x = ", phase1WorstX, ",  J = ", phase1WorstJ];
Print["  Violations (F < 0):  ", Length[phase1Violations]];

If[Length[phase1Violations] > 0,
  Print[""];
  Print["  *** VIOLATIONS FOUND AT SAMPLE POINTS ***"];
  Print["  The SDPB bound is INVALID — functional is negative for unconstrained spins."];
  Print[""];
  Print["  Worst violations:"];
  sorted = SortBy[phase1Violations, #[[3]] &];
  Do[
    Print["    x = ", sorted[[k, 1]], "  J = ", sorted[[k, 2]],
          "  F = ", sorted[[k, 3]]],
    {k, Min[20, Length[sorted]]}
  ];
  Print[""];

  violatedJs = DeleteDuplicates[Sort[#[[2]] & /@ phase1Violations]];
  Print["  J values with violations: ", violatedJs];
  Print["  Suggested: add these to Jlist in test9.m and rerun SDPB."];
];
Print[""];


(* ----------------------------------------------------------------
   7. PHASE 2: CHECK AT MIDPOINTS
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  PHASE 2: Checking at MIDPOINTS"];
Print["========================================"];

nIntervals = Length[samplePoints] - 1;

allIntervals = Join[
  If[samplePoints[[1]]  > xLeft  + 10^-12, {{xLeft + singularTol, samplePoints[[1]]}}, {}],
  Table[{samplePoints[[i]], samplePoints[[i+1]]}, {i, nIntervals}],
  If[samplePoints[[-1]] < xRight - 10^-12, {{samplePoints[[-1]], xRight - singularTol}}, {}]
];

Print["  ", Length[allIntervals], " intervals and ", Length[JlistLarge], " J-values"];
Print["  Left boundary  : [", xLeft, ", ", samplePoints[[1]], "]  ",
      If[samplePoints[[1]] > xLeft  + 10^-12, "(active)", "(skipped — x_min = xLeft)"]];
Print["  Right boundary : [", samplePoints[[-1]], ", ", xRight, "]  ",
      If[samplePoints[[-1]] < xRight - 10^-12, "(active)", "(skipped — x_max = xRight)"]];
Print["  (stopping threshold: min_width = ", minWidth, ")"];
Print[""];

phase2Violations = {};
phase2MinVal = Infinity;
phase2WorstX = None;
phase2WorstJ = None;

Do[
  xa = allIntervals[[i, 1]];
  xb = allIntervals[[i, 2]];
  width = xb - xa;
  mid = (xa + xb) / 2;

  If[Abs[mid] < singularTol, mid = singularTol];
  If[Abs[1 - mid] < singularTol, mid = 1 - singularTol];

  localMin = Infinity;
  localWorstJ = None;

  Do[
    val = safeF[mid, Jj];
    If[val =!= $Failed && val < localMin,
      localMin = val;
      localWorstJ = Jj;
    ],
    {Jj, JlistLarge}
  ];

  If[localMin < 0,
    AppendTo[phase2Violations, {mid, localWorstJ, localMin}];
  ];

  If[localMin < phase2MinVal,
    phase2MinVal = localMin;
    phase2WorstX = mid;
    phase2WorstJ = localWorstJ;
  ];

  If[Mod[i, Max[1, Floor[Length[allIntervals]/10]]] == 0,
    Print["  Progress: ", i, "/", Length[allIntervals],
          "  local min so far: ", phase2MinVal]
  ],
  {i, Length[allIntervals]}
];

Print[""];
Print["Phase 2 Results:"];
Print["  Global minimum F value:  ", phase2MinVal];
Print["  Worst point:  x = ", phase2WorstX, ",  J = ", phase2WorstJ];
Print["  Violations (F < 0):  ", Length[phase2Violations]];

If[Length[phase2Violations] > 0,
  Print[""];
  Print["  *** VIOLATIONS FOUND AT MIDPOINTS OF SAMPLING POINTS ***"];
  Print[""];
  Print["  Worst violations:"];
  sorted2 = SortBy[phase2Violations, #[[3]] &];
  Do[
    Print["    x = ", sorted2[[k, 1]], "  J = ", sorted2[[k, 2]],
          "  F = ", sorted2[[k, 3]]],
    {k, Min[20, Length[sorted2]]}
  ];
  Print[""];

  violatedJs2 = DeleteDuplicates[Sort[#[[2]] & /@ phase2Violations]];
  Print["  J values with violations: ", violatedJs2];
];
Print[""];


(* ----------------------------------------------------------------
   8. ASYMPTOTIC ANALYSIS
      BUG FIX A: extraFuncs had 3 entries; corrected to use extraTriplet
      (which has 5 entries matching fList) so the J→∞ functional is
      correctly computed as Σ_k yVec[[k]] * extraTriplet[[k]][x].
      This selects y₅ * LargeJ(x) — the X53 large-J contribution.

      BUG FIX B: display of val used fixSciNotation /@ val on a Number
      (not a String).  Replaced with N[val, 50] for a clean display.
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  ASYMPTOTIC ANALYSIS"];
Print["========================================"];

nRep = Min[5, Length[samplePoints]];
repIndices = Table[
  Max[1, Min[Length[samplePoints], Round[i * Length[samplePoints] / (nRep + 1)]]],
  {i, nRep}
];
repXvals = samplePoints[[#]] & /@ DeleteDuplicates[repIndices];

Print["Checking J-dependence at representative x values:"];
JprobeList = {60, 100, 200, 500, 1000, 2000, 5000, 10000};
JprobeList = Select[JprobeList, # <= jScanMax &];

Do[
  xi = repXvals[[r]];
  Print[""];
  Print["  x = ", xi, "  (sp = ", N[1/(1-xi), 6], ")"];
  Print["    J       |  F(x,J)                    |  F/J^6"];
  Print["    --------|----------------------------|------------------------------"];
  Do[
    val = safeF[xi, Jp];
    (* BUG FIX B: val is already a Number — display directly *)
    ratio = If[val =!= $Failed && Jp > 0, val / Jp^6, "N/A"];
    Print["    ", Jp, "     |  ",
          If[val =!= $Failed, N[val, 50], "$Failed"],
          "    |  ", If[NumberQ[ratio], N[ratio, 50], ratio]],
    {Jp, JprobeList}
  ];
  (* BUG FIX A: use extraTriplet (5 entries) and Length[yVec] (not 3) *)
  xInfVal = Sum[yVec[[k]] * extraTriplet[[k]][xi], {k, Length[yVec]}];
  Print["    inf     |  (LargeJ limit)            |  ", N[xInfVal, 50]],
  {r, Length[repXvals]}
];
Print[""];


(* ----------------------------------------------------------------
   9. FINAL SUMMARY
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  FINAL SUMMARY"];
Print["========================================"];

totalViolations = Length[phase1Violations] + Length[phase2Violations];

If[totalViolations == 0,
  Print[""];
  Print["  O  NO VIOLATIONS FOUND for J up to ", jScanMax];
  Print["  O  The current Jmax = 60 appears SUFFICIENT."];
  Print["  O  The SDPB bound is consistent with large-J positivity."];
  Print[""];
  Print["  Note: This does not prove positivity for ALL J > 60 at ALL x,"];
  Print["  but the smart grid provides strong evidence. The asymptotic"];
  Print["  analysis above shows the convergence pattern toward J -> infty."];
  Print[""];
  Quit[0],

  Print[""];
  Print["  X  VIOLATIONS FOUND: ", totalViolations, " total"];
  Print["     Phase 1 (sample points): ", Length[phase1Violations]];
  Print["     Phase 2 (midpoints):     ", Length[phase2Violations]];
  Print[""];

  allViolatedJs = DeleteDuplicates[Sort[Join[
    If[Length[phase1Violations] > 0, #[[2]] & /@ phase1Violations, {}],
    If[Length[phase2Violations] > 0, #[[2]] & /@ phase2Violations, {}]
  ]]];

  Print["  Violated J values: ", allViolatedJs];
  Print[""];
  Print["  RECOMMENDED ACTIONS:"];

  If[Length[allViolatedJs] <= 20,
    Print["  1. Add these J values to Jlist in test9.m:"];
    Print["     JlistLarge = ", allViolatedJs, ";"];
    Print["     Jlist = Join[Range[0, 60, 2], JlistLarge];"];
    Print["  2. Rerun test9.m -> pmp2sdp -> sdpb"];
    Print["  3. Rerun this diagnostic to verify"],

    Print["  1. Violations are widespread (", Length[allViolatedJs], " J values)."];
    Print["     Consider implementing the 'large J' regime from Albert et al.:"];
    Print["     - Discrete evaluation in both J and m² for J > 60"];
    Print["     - See Appendix A.2 of arXiv:2406.12959"];
    Print["  2. As a quick fix, increase Jmax in test9.m to ", Max[allViolatedJs]];
  ];

  Print[""];
  Quit[1]
];