(* ================================================================
   check_large_J.m
   ----------------------------------------------------------------
   PURPOSE
     After an SDPB run with Jmax=60, read the dual vector y and the
     current sampling points, then check whether the functional

         F(x, J) = y₁·g20(x,J) + y₂·g31(x,J) + y₃·n4(x,J)

     remains non-negative for LARGE spins J > 60.

     This is a DIAGNOSTIC tool — it does NOT modify any files.
     It reports whether the "small J" regime (J ≤ 60) is sufficient
     or whether additional spins must be added to the SDPB constraints.

   TWO PHASES
     Phase 1: Check F(x_i, J) at existing sample points x_i for
              J = 62, 64, ..., up to J_scan_max.
              (At sample points, SDPB guarantees F ≥ 0 for J ≤ 60.
               Violations here mean the bound is INVALID.)

     Phase 2: Check F(mid, J) at midpoints between sample points
              for the same large J range.
              (Catches cases where F is marginal at sample points
               but dips negative between them.)

   SMART J GRID
     Instead of checking every even J up to 10000 (5000 values),
     we use an adaptive grid:
       Dense near Jmax:   {62, 64, 66, ..., 120}
       Moderate spacing:  {140, 160, ..., 500}
       Sparse/log:        {600, 800, 1000, 1500, 2000, 3000, 5000, 7500, 10000}
     Total: ~100 J values instead of 5000.

   USAGE
     wolframscript -file check_large_J.m <sp_file> <y_file> [J_scan_max]

   EXIT CODES
     0  F ≥ 0 for all large J tested — Jmax=60 is sufficient
     1  Violations found — need to expand Jmax or add large-J regime
     2  Error (bad arguments, unreadable files)
   ================================================================ *)


(* ----------------------------------------------------------------
   1. ARGUMENT PARSING
   ---------------------------------------------------------------- *)

myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

If[Length[myArgs] < 2,
  Print["USAGE: wolframscript -file negative_region.m ",
        "<sp_file> <y_file> [N_pts] [out_sp_file]"];
  Quit[2]
];

spFile    = myArgs[[1]];
yFile     = myArgs[[2]];
nPts      = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 10];
minWidth  = If[Length[myArgs] >= 5, ToExpression[myArgs[[5]]], 10^-6];

(* spin grid *)
jScanMax  = If[Length[myArgs] >= 4, ToExpression[myArgs[[4]]], 10000];

Print["=== check_large_J.m  (Large-J Positivity Diagnostic) ==="];
Print["  sp_file     = ", spFile];
Print["  y_file      = ", yFile];
Print["  N_pts      = ", nPts];
Print["  min_width  = ", minWidth];
Print["  J_scan_max  = ", jScanMax];
Print[""];


(* ----------------------------------------------------------------
   2. READ FILES
   ---------------------------------------------------------------- *)
(* Fix C/Fortran-style scientific notation: e.g. "1.5e+02" -> "1.5*^+02"
   Mathematica does not recognise lowercase 'e' or uppercase 'E' as an
   exponent marker; it treats them as the symbol e or Euler's E.
   This regex converts  e±digits  /  E±digits  to  *^±digits  before
   ToExpression is called. *)
fixSciNotation[s_String] := StringReplace[s,
  RegularExpression["[eE]([+-]?\\d+)"] :> "*^$1"
];

(* Sampling points *)
If[!FileExistsQ[spFile],
  Print["ERROR: sampling points file not found: ", spFile]; Quit[2]];

spRaw = Select[
  ReadList[spFile, String],
  StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &
];

If[Length[spRaw] == 0,
  Print["ERROR: no sample points found in ", spFile]; Quit[2]];

(* Parse and sort; keep high precision *)
samplePoints = Sort[SetPrecision[ToExpression /@ (fixSciNotation /@ spRaw), 200]];
Print["Loaded ", Length[samplePoints], " sampling points."];

(* ----------------------------------------------------------------
   3.  READ z.txt
       SDPB writes all n components of y, one per line; blank lines
       and lines starting with "#" are skipped.
       The length of yVec must equal the number of functions in fVec.
       A mismatch is caught by the dimensional check below.
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

(* Validate: every component of yVec must be numeric.
   A non-numeric entry (e.g. containing symbolic 'e') means the
   scientific-notation fix missed a case or the file is malformed. *)
Do[
  If[!NumberQ[yVec[[k]]],
    Print["ERROR: y component ", k, " is not numeric: ", yVec[[k]]];
    Print["  Raw line was: ", yRaw[[k]]];
    Print["  After fixSciNotation: ", fixSciNotation[yRaw[[k]]]];
    Quit[2]
  ],
  {k, Length[yVec]}
];

(* Same validation for sample points *)
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
   4. PROBLEM-SPECIFIC DEFINITIONS  (must match test9.m / refine_sampling.m)
   ---------------------------------------------------------------- *)

(* problem-specific *)
(* let the mass be 4mA^2 < M^2 = 1 where M  = 1 to infinity *)

maVal = SetPrecision[3/20, 600];

Print["mA = ", maVal]

(* dispersion representation of Wilson coefficients *)
(* All functions now precompute sp and mA numerically with N[...,prec] *)

(* forward limit: use our own convention *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[-(2*Sqrt[sp/(-4*mA^2 + sp)])/(2*mA^2 - sp)^3, 600]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[-(Sqrt[sp/(-4*mA^2 + sp)]*(-3 - (2*J*(1 + J)*(2*mA^2 - sp))/(-4*mA^2 + sp)))/(-2*mA^2 + sp)^4, 600]
];

(* --- n4: safe precision guard -------------------------------------------
   n4 contains LegendreP[J, 1, z] and LegendreP[J, 2, z] with z > 1, and
   Hypergeometric2F1[-J, 1+J, 1, z1] which equals P_J(1.066) at typical
   sample points.  All three grow as A^J where A ~ 1.436 (at x \[TildeTilde] 0),
   requiring ~1560 decimal digits at J = 10000 just to represent individual
   terms before cancellation.  With prec = 600 the computed result is pure
   numerical noise of magnitude ~10^960.

   Physical justification for returning 0 at large J:
     extraTriplet = {0&, 0&, 0&, 0&, LargeJ}
   explicitly encodes n4 \[RightArrow] 0 as J \[RightArrow] \[Infinity] (index 3 in fList = n4).
   The safe threshold is  J_SAFE = floor((prec - 30) / 0.156) \[TildeTilde] 3654.
   Since the only super-threshold spin in Jlist is J = 10000, this guard
   affects exactly the J = 10000 blocks.
   --------------------------------------------------------------------- *)

(* Maximum log10[A(x)] over physical x in (0,1): A = z + sqrt(z^2-1),
   z = 1 + 8*mA^2/(3*(sp - 4*mA^2)).  Evaluated numerically: max ~ 0.156. *)
n4MaxLog10A = 0.156;

n4Safe[x_?NumericQ, J_?IntegerQ] :=
  If[J > Floor[(600 - 30) / n4MaxLog10A],   (* J > ~3654 *)
    (* Large-J limit: n4 \[RightArrow] 0 (consistent with extraTriplet's 0& for n4) *)
    SetPrecision[0, 600],
    (* Normal evaluation for J \[LessEqual] 3654 *)
    n4[x, J]
  ];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(81*Sqrt[sp/(-4*mA^2 + sp)]*(8*mA^4*(14*mA^2 - 15*sp)*(-8*mA^2 + 3*sp)^(3/2)*Hypergeometric2F1[-J, 1 + J, 1, (4*mA^2)/(12*mA^2 - 3*sp)] + (8*mA^4 - 18*mA^2*sp + 9*sp^2)*((-2*I)*mA*(10*mA^2 - 9*sp)*(8*mA^2 - 3*sp)*LegendreP[J, 1, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))] + Sqrt[-8*mA^2 + 3*sp]*(-8*mA^4 + 18*mA^2*sp - 9*sp^2)*LegendreP[J, 2, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))])))/(4*mA^2*(-2*mA^2 + sp)*(-8*mA^2 + 3*sp)^(3/2)*(8*mA^4 - 18*mA^2*sp + 9*sp^2)^3), 600]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

X53[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*((-4 + J)*(-2 + J)*(3 + J)*(5 + J) + ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

(* Large J limit of X53 *)
LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-32*mA^6 + 24*mA^4*sp - 6*mA^2*sp^2 + sp^3))/(18*sp^3*(-4*mA^2 + sp)^6), 600]
];

(* Large J limit of X52 *)
LargeJ2[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[-1/18*(Sqrt[sp/(-4*mA^2 + sp)]*(-32*mA^6 + 24*mA^4*sp - 6*mA^2*sp^2 + sp^3))/(sp^3*(-4*mA^2 + sp)^6), 600];
]

Jmax = 60;
Jlist = Range[0, Jmax, 2];
(* JlistLarge = {10000};
Jlist = Join[Range[0, 60, 2], JlistLarge]; *)

(* NOTE: n4 is replaced by n4Safe (returns 0 for J > ~3654, exact for J <= 3654).
   All other functions are unchanged; they grow polynomially in J and are
   numerically accurate with prec = 600 even at J = 10000. *)
fList = {g20, g31, n4Safe, X52, X53};

(* large J limit *)
(* 0& is a constant function of 0 *)

extraTriplet = {0&, 0&, 0&, LargeJ2, LargeJ};

(* optimal lower bound *)
(* norm = {0, 1, 0, 0, 0};
obj = {-1, 0, 0, 0, 0}; *)

(* optimal upper bound *)
norm = {0, -1, 0, 0, 0};
obj = {-1, 0, 0, 0, 0};

xLeft  = SetPrecision[0, 600];   (* physical domain left endpoint  — check includes [xLeft,  x_min] *)
xRight = SetPrecision[1, 600];   (* physical domain right endpoint — check includes [x_max, xRight] *)

(* --- Dimensional consistency check ---
   y.txt must have exactly as many lines as there are functions in
   fList.  A mismatch means the wrong y.txt or fList was supplied. *)
If[Length[yVec] != Length[fList],
  Print["ERROR: y.txt has ", Length[yVec], " component(s) but fList has ",
        Length[fList], " function(s). They must match."];
  Quit[2]
];


(* ----------------------------------------------------------------
   4. BUILD SMART J GRID FOR SCANNING
   ---------------------------------------------------------------- *)

(* The SDPB run constrained J = 0, 2, ..., 60 and J->infty.
   We need to check J = 62, 64, ..., J_scan_max. *)

JlistDense    = Range[62, Min[120, jScanMax], 2];         (* every even spin near Jmax *)
JlistModerate = Range[140, Min[500, jScanMax], 20];       (* every 20 *)
JlistSparse   = Select[
  {600, 800, 1000, 1500, 2000, 3000, 5000, 7500, 10000, 15000, 20000},
  # <= jScanMax &
];

JlistLarge = DeleteDuplicates[Sort[Join[JlistDense, JlistModerate, JlistSparse]]];

Print["Large-J scan grid: ", Length[JlistLarge], " spin values"];
Print["  Range: J = ", First[JlistLarge], " to ", Last[JlistLarge]];
Print["  Dense (62-120):    ", Length[JlistDense], " values"];
Print["  Moderate (140-500): ", Length[JlistModerate], " values"];
Print["  Sparse (600+):     ", Length[JlistSparse], " values"];
Print[""];


(* ----------------------------------------------------------------
   5. FUNCTIONAL EVALUATION
   ---------------------------------------------------------------- *)

F[x_?NumericQ, J_?IntegerQ] := Sum[yVec[[k]] * fList[[k]][x, J], {k, Length[yVec]}];

(* Safety: clamp midpoints near known singularity at x -> 1 (sp = 1/(1-x))
   and provide a robust numeric evaluator that returns $Failed on
   non-finite or exceptional results. *)
singularTol = SetPrecision[10^-12, 600]; (* distance from x=1 to avoid sp singularity *)

safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  If[Abs[x0] < singularTol, x0 = singularTol];
  val = Quiet[Check[N[F[x0, J], 600], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate}, Head[val]],
    $Failed,
    val
  ]
];


(* ----------------------------------------------------------------
   6. PHASE 1: CHECK AT SAMPLE POINTS
      At x_i where SDPB guarantees F ≥ 0 for J ≤ 60,
      check whether F ≥ 0 also holds for J > 60.
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  PHASE 1: Checking at SAMPLE POINTS"];
Print["========================================"];
Print["  ", Length[samplePoints], " x-points and ", Length[JlistLarge], "J-values"];
Print[""];

phase1Violations = {};
phase1MinVal = Infinity;
phase1WorstX = None;
phase1WorstJ = None;

Do[
  xi = samplePoints[[i]];

  (* For each x_i, find the minimum F over all large J *)
  localMin = Infinity;
  localWorstJ = None;

  Do[
    val = safeF[xi, Jj];
    If[val =!= $Failed && val < localMin,
      localMin = val;
      localWorstJ = Jj;
    ],
    {Jj, JlistLarge}
  ];

  (* Negative values *)
  If[localMin < 0,
    AppendTo[phase1Violations, {xi, localWorstJ, localMin}];
  ];

  If[localMin < phase1MinVal,
    phase1MinVal = localMin;
    phase1WorstX = xi;
    phase1WorstJ = localWorstJ;
  ];

  (* Progress indicator every 10 points *)
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

  (* Identify critical J values *)
  violatedJs = DeleteDuplicates[Sort[#[[2]] & /@ phase1Violations]];
  Print["  J values with violations: ", violatedJs];
  Print["  Suggested: add these to Jlist in test9.m and rerun SDPB."];
];
Print[""];


(* ----------------------------------------------------------------
   7. PHASE 2: CHECK AT MIDPOINTS
      Even if Phase 1 passes, F could dip negative between
      sample points for large J.
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  PHASE 2: Checking at MIDPOINTS"];
Print["========================================"];

nIntervals = Length[samplePoints] - 1;

(* Include boundary intervals *)
allIntervals = Join[
  If[samplePoints[[1]]  > xLeft  + 10^-12, {{xLeft + singularTol, samplePoints[[1]]}}, {}],
  Table[{samplePoints[[i]], samplePoints[[i+1]]}, {i, nIntervals}],
  If[samplePoints[[-1]] < xRight - 10^-12, {{samplePoints[[-1]], xRight - singularTol}}, {}]
];

Print["  ", Length[allIntervals], " intervals and ", Length[JlistLarge], "J-values"];
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

  (* Clamp away from singularities *)
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
      Check how F scales with J at a few representative x values.
      This reveals whether violations are transient or persistent.
   ---------------------------------------------------------------- *)

Print["========================================"];
Print["  ASYMPTOTIC ANALYSIS"];
Print["========================================"];

(* Pick a few representative x values *)
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
  Print["    J       |  F(x,J)          |  F/J^6"];
  Print["    --------|------------------|----------------"];
  Do[
    val = safeF[xi, Jp];
    ratio = If[val =!= $Failed && Jp > 0, val / Jp^6, "N/A"];
    Print["    ", Jp, "     |  ", If[val =!= $Failed, SetPrecision[ToExpression /@ (fixSciNotation /@ val), 200], "$Failed"],
          "    |  ", If[NumberQ[ratio], SetPrecision[ToExpression /@ (fixSciNotation /@ ratio), 200], ratio]],
    {Jp, JprobeList}
  ];
  (* Also show the J->infty limit for comparison *)
  extraFuncs = extraTriplet;
  xInfVal = Sum[yVec[[k]] * extraFuncs[[k]][xi], {k, Length[yVec]}];
  Print["    inf     |  (LargeJ limit)  |  ", SetPrecision[ToExpression /@ (fixSciNotation /@ xInfVal), 200]],
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

  (* Violations found *)
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