(* ================================================================
   refine_sampling.m  (precision-propagation corrected version)
   ----------------------------------------------------------------
   PURPOSE
     After each SDPB run, read the current sampling points and the
     solver's y.txt, reconstruct the functional

         F(x, J) = f₁(x,J)·y₁ + f₂(x,J)·y₂ + f₃(x,J)·y₃ + ...

     for EVERY J in Jlist (discrete spin, constrained exactly), and
     identify x-intervals where F is negative at the midpoint for
     any J.

   PRECISION FIXES (2026-05)
     Applied consistently with test9.m and largeSpinTest.m:
     1. maVal: N[3/20, 650] instead of SetPrecision[0.150, 600].
     2. n4PrecMin: adaptive precision for n4 based on J and x.
     3. n4: evaluates sp, mA at prec_local digits; Re[…] explicit.
     4. All other functions: N[3/20, 650] for mA.

     For J ≤ 60 (normal Jlist): prec_local = 650, zero overhead.
   ================================================================ *)


(* ----------------------------------------------------------------
   1.  ARGUMENT PARSING
   ---------------------------------------------------------------- *)

myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

If[Length[myArgs] < 2,
  Print["USAGE: wolframscript -file refine_sampling.m ",
        "<sp_file> <y_file> [N_pts] [out_sp_file]"];
  Quit[2]
];

spFile    = myArgs[[1]];
yFile     = myArgs[[2]];
nPts      = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 10];
minWidth  = If[Length[myArgs] >= 4, ToExpression[myArgs[[4]]], 10^-6];
outFile   = If[Length[myArgs] >= 5, myArgs[[5]], spFile];

Print["=== refine_sampling.m ==="];
Print["  sp_file    = ", spFile];
Print["  y_file     = ", yFile];
Print["  N_pts      = ", nPts];
Print["  min_width  = ", minWidth];
Print["  out_file   = ", outFile];
Print[""];


(* ----------------------------------------------------------------
   2.  READ SAMPLING POINTS
   ---------------------------------------------------------------- *)

If[!FileExistsQ[spFile],
  Print["ERROR: sampling points file not found: ", spFile]; Quit[2]];

spRaw = Select[
  ReadList[spFile, String],
  StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &
];

If[Length[spRaw] == 0,
  Print["ERROR: no sample points found in ", spFile]; Quit[2]];

fixSciNotation[s_String] := StringReplace[s,
  RegularExpression["[eE]([+-]?\\d+)"] :> "*^$1"
];

samplePoints = Sort[SetPrecision[ToExpression /@ (fixSciNotation /@ spRaw), 200]];

Print["Loaded ", Length[samplePoints], " sampling points:"];
Print["  ", samplePoints];
Print[""];


(* ----------------------------------------------------------------
   3.  READ y.txt
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


(* ================================================================
   PROBLEM-SPECIFIC SECTION — must match test9.m
   ================================================================ *)

(* FIX 1: exact rational for mA — avoids machine-float precision leakage *)
maVal = N[3/20, 650];

Print["mA = ", maVal];

(* FIX 2: n4PrecMin — adaptive working precision for n4.
   For J ≤ ~3782: returns 650 (no overhead).
   For J > 3782:  returns J * 0.1568 + 60 (precision required to resolve
   the catastrophic cancellation of the A(x)^J terms in n4). *)
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
   FIX 2+3+4: adaptive precision, exact-rational mA, explicit Re[…] *)
n4[x_?NumericQ, J_?IntegerQ] := Module[{prec, xp, sp, mA, result},
  prec = n4PrecMin[J, x];
  xp   = SetPrecision[x, prec];    (* promote x to avoid low-precision bottleneck *)
  sp   = N[1/(1 - xp), prec];
  mA   = N[maVal, prec];            (* exact rational, independent of global maVal *)
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

extraTriplet = {0&, 0&, 0&, 0&, LargeJ};

norm = {0, -1, 0, 0, 0};
obj  = {-1, 0, 0, 0, 0};

xLeft  = SetPrecision[0, 650];
xRight = SetPrecision[1, 650];

(* ================================================================
   END OF PROBLEM-SPECIFIC SECTION
   ================================================================ *)


If[Length[yVec] != Length[fList],
  Print["ERROR: y.txt has ", Length[yVec], " component(s) but fList has ",
        Length[fList], " function(s). They must match."];
  Quit[2]
];

If[Length[extraTriplet] != Length[fList],
  Print["ERROR: extraTriplet has ", Length[extraTriplet],
        " entries but fList has ", Length[fList], ". They must match."];
  Quit[2]
];

F[x_?NumericQ, J_?IntegerQ] := Sum[yVec[[k]] * fList[[k]][x, J], {k, Length[yVec]}];
X[x_?NumericQ] := Sum[yVec[[k]] * extraTriplet[[k]][x], {k, Length[yVec]}];

singularTol = 10^-12;

safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  val = Quiet[Check[N[F[x0, J], 200], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];

safeX[x_?NumericQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  val = Quiet[Check[N[X[x0], 200], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];


(* ----------------------------------------------------------------
   5.  MIDPOINT CHECK OVER CONSECUTIVE INTERVALS
   ---------------------------------------------------------------- *)

nIntervals    = Length[samplePoints] - 1;
interiorPairs = Table[
  {samplePoints[[i]], samplePoints[[i + 1]]},
  {i, nIntervals}
];
allIntervals = Join[
  If[samplePoints[[1]]  > xLeft  + 10^-12, {{xLeft,  samplePoints[[1]]}},  {}],
  interiorPairs,
  If[samplePoints[[-1]] < xRight - 10^-12, {{samplePoints[[-1]], xRight}}, {}]
];
nAllIntervals    = Length[allIntervals];
flaggedIntervals = {};
tinyNegative     = {};
midpointValues   = {};

Print["Checking ", nAllIntervals, " interval(s) on [", xLeft, ", ", xRight,
      "]  \[Times]  ", Length[Jlist], " J-value(s):"];
Print["  Interior pairs : ", Length[interiorPairs]];
Print["  Left boundary  : [", xLeft, ", ", samplePoints[[1]], "]  ",
      If[samplePoints[[1]] > xLeft  + 10^-12, "(active)", "(skipped — x_min = xLeft)"]];
Print["  Right boundary : [", samplePoints[[-1]], ", ", xRight, "]  ",
      If[samplePoints[[-1]] < xRight - 10^-12, "(active)", "(skipped — x_max = xRight)"]];
Print["  (stopping threshold: min_width = ", minWidth, ")"];
Print[""];

Do[
  xa    = allIntervals[[i, 1]];
  xb    = allIntervals[[i, 2]];
  width = xb - xa;
  mid   = (xa + xb) / 2;

  fmByJ = Table[safeF[mid, Jlist[[j]]], {j, Length[Jlist]}];
  failedCount = Count[fmByJ, $Failed];
  If[failedCount > 0,
    Print["  WARNING: ", failedCount,
          " non-finite F(mid,J) evaluation(s) at mid=", mid,
          " — proceeding with available J values (if any)."]
  ];

  xVal = safeX[mid];
  If[xVal === $Failed,
    Print["  NOTE: X(mid) evaluation failed at mid=", mid]
  ];

  validPairs = Select[Transpose[{Jlist, fmByJ}], Last[#] =!= $Failed &];

  If[validPairs == {} && xVal === $Failed,
    Print["  ERROR: F(mid,J) failed for all J and X(mid) failed at mid=", mid, "; flagging."];
    AppendTo[flaggedIntervals, {xa, xb}];
    AppendTo[midpointValues, {mid, $Failed, "all_failed"}];
    Continue[];
  ];

  If[validPairs == {},
    fmJ = $Failed,
    vals = Last /@ validPairs;
    js   = First /@ validPairs;
    fmJ  = Min[vals];
    worstJIndex     = First[Ordering[vals, 1]];
    worstJcandidate = js[[worstJIndex]];
  ];

  Which[
    fmJ === $Failed && xVal =!= $Failed,
      fm = xVal; worstSource = "X"; worstJ = None,

    fmJ =!= $Failed && xVal === $Failed,
      fm = fmJ; worstSource = "J"; worstJ = worstJcandidate,

    fmJ =!= $Failed && xVal =!= $Failed,
      If[xVal <= fmJ, fm = xVal; worstSource = "X"; worstJ = None,
         fm = fmJ; worstSource = "J"; worstJ = worstJcandidate
      ]
  ];

  AppendTo[midpointValues, {mid, fm, worstSource}];

  Which[
    fm >= 0,
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  min value=", fm, "  ok (source=", worstSource, ")"],

    fm < 0 && width < minWidth,
      AppendTo[tinyNegative, {xa, xb}];
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  min value=", fm, "  (source=", worstSource, ")",
            "  negative but w < ", minWidth, " — below threshold, skipping"],

    fm < 0 && width >= minWidth,
      AppendTo[flaggedIntervals, {xa, xb}];
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  min value=", fm, "  (source=", worstSource,
            If[worstSource=="J", ", worst J=", ""],
            If[worstSource=="J", worstJ, ""], ")",
            "  *** NEGATIVE, will refine ***"]
  ],
  {i, nAllIntervals}
];
Print[""];


(* ----------------------------------------------------------------
   6.  CONVERGENCE CHECK
   ---------------------------------------------------------------- *)

If[Length[flaggedIntervals] == 0,
  If[Length[tinyNegative] == 0,
    Print["CONVERGED: F(midpoint) >= 0 for all ", nAllIntervals, " intervals."],
    Print["CONVERGED: ", Length[tinyNegative],
          " interval(s) have F(mid) < 0 but all widths < ", minWidth, "."];
    Print["  Tiny negative intervals (not refined):"];
    Do[Print["    [", r[[1]], ", ", r[[2]], "]  width = ", r[[2]] - r[[1]]],
       {r, tinyNegative}]
  ];
  Print["No new sampling points needed. ", outFile, " is unchanged."];
  Quit[0]
];

Print[Length[flaggedIntervals], " refinable interval(s). Generating new points..."];
If[Length[tinyNegative] > 0,
  Print["  (", Length[tinyNegative],
        " additional interval(s) negative but below threshold — not refined)"]
];
Print[""];


(* ----------------------------------------------------------------
   7.  GENERATE NEW SAMPLING POINTS
   ---------------------------------------------------------------- *)

newPoints = Flatten[
  Table[
    Module[{xa, xb, xStar, s},
      xa    = interval[[1]];
      xb    = interval[[2]];
      If[Abs[xa - xLeft] < singularTol, xa = xLeft + singularTol];
      If[Abs[xb - xRight] < singularTol, xb = xRight - singularTol];
      xStar = (xa + xb) / 2;
      s     = (xb - xa) / nPts;
      Print["  Flagged interval [", interval[[1]], ", ", interval[[2]], "] (clamped to [", xa, ", ", xb, "]):"];
      Print["    x* = ", xStar, "  s = ", s];
      Table[SetPrecision[xStar + (k - nPts/2) * s, 200], {k, 0, nPts}]
    ],
    {interval, flaggedIntervals}
  ]
];

Print[""];
Print["New candidate points (before dedup): ", Length[newPoints]];

dedupTol = 10^-10;


(* ----------------------------------------------------------------
   8.  MERGE AND DEDUPLICATE
   ---------------------------------------------------------------- *)

allPoints = Sort[Join[samplePoints, newPoints]];

dedupPoints = Fold[
  Function[{kept, pt},
    If[Abs[pt - Last[kept]] < dedupTol, kept, Append[kept, pt]]
  ],
  {First[allPoints]},
  Rest[allPoints]
];

dedupPoints = Select[dedupPoints,
  Function[x, Abs[x - xLeft] >= singularTol && Abs[x - xRight] >= singularTol]
];

Print["Original points          : ", Length[samplePoints]];
Print["New refined points       : ", Length[newPoints]];
Print["Combined before dedup    : ", Length[allPoints]];
Print["Combined after  dedup    : ", Length[dedupPoints]];
Print[""];


(* ----------------------------------------------------------------
   9.  WRITE MERGED SAMPLING POINTS TO FILE
   ---------------------------------------------------------------- *)

formatPlainDecimal[x_?NumericQ] := Module[
  {ax, sign, digits, dpos, ndig = 50, dstr, result},
  If[x == 0, Return["0." <> StringJoin[Table["0", {ndig}]]]];
  sign = If[Negative[x], "-", ""];
  ax   = SetPrecision[Abs[x], ndig];
  {digits, dpos} = RealDigits[ax, 10, ndig];
  dstr = StringJoin[ToString /@ Replace[digits, n_ /; n < 0 :> 0, {1}]];
  result = Which[
    dpos >= ndig,
      dstr <> StringJoin[Table["0", {dpos - ndig}]] <> ".0",
    dpos > 0,
      StringTake[dstr, dpos] <> "." <> StringDrop[dstr, dpos],
    dpos == 0,
      "0." <> dstr,
    dpos < 0,
      "0." <> StringJoin[Table["0", {-dpos}]] <> dstr
  ];
  sign <> result
];

outLines = StringRiffle[
  formatPlainDecimal /@ dedupPoints,
  "\n"
];

Export[outFile, outLines, "Text"];

Print["Wrote ", Length[dedupPoints], " sampling points to: ", outFile,
      "  (", Length[samplePoints], " original + ",
      Length[dedupPoints] - Length[samplePoints], " new refined)"];
Print[""];
Print["New sampling points:"];
Print["  ", dedupPoints];
Print[""];
Print["STATUS: Refinement needed. Re-run pmp2sdp + sdpb with updated ", outFile, "."];

Quit[1];