(* ================================================================
   negative_region.m
   ----------------------------------------------------------------
   PURPOSE
     After each SDPB run, read the current sampling points and the
     solver's y.txt, reconstruct the functional

         F(x, J) = f₁(x,J)·y₁ + f₂(x,J)·y₂ + f₃(x,J)·y₃

     for EVERY J in Jlist (discrete spin, constrained exactly), and
     identify x-intervals where F is negative at the midpoint for
     any J.

   WHY MIDPOINTS ONLY (not a fine grid scan)
     SDPB already guarantees F(xᵢ, J) >= 0 at every sampling point
     xᵢ and for every J in Jlist.  Therefore we do NOT need to scan
     all x.  We probe each open interval (xᵢ, xᵢ₊₁) by evaluating
     F(mid, J) for ALL J in Jlist.  If the minimum over J is negative
     the interval is flagged.  This is faster and more principled
     than a fine-grid sweep.

   J IS NOT SAMPLED — IT IS EXACT
     x is continuous and is discretised adaptively.
     J is discrete; every J in Jlist is checked at each midpoint.
     This is x-independent and is
     checked once before the interval loop; x-refinement cannot
     fix a violation of this constraint.

   NEW POINTS (from observation.md)
     For each flagged interval [xa, xb]:
       x*  = (xa + xb) / 2          (midpoint = centre of interval)
       s   = (xb - xa) / N_pts      (step size)
       new = { x* + (k - N_pts/2)*s : k = 0, 1, ..., N_pts }
           = { xa, xa+s, ..., x*, ..., xb-s, xb }
     This places N_pts+1 points that span [xa, xb] exactly.
     xa and xb are already sample points; after deduplication they
     are dropped, leaving N_pts-1 new interior points.

   STOPPING CRITERION
     If F(midpoint) < 0 but the interval width  xb - xa < min_width
     (default 10^-6), the interval is NOT flagged for refinement.
     The functional is negative there, but the region is so narrow
     that adding more points would not improve the bound meaningfully.
     Convergence is declared when every interval with F(mid) < 0 has
     width below min_width — i.e. no refinable negative intervals remain.

   USAGE
     wolframscript -file negative_region.m <sp_file> <y_file> \
                   [N_pts] [min_width] [out_sp_file]

     sp_file      path to current sampling_points.txt (one x per line)
     y_file       path to SDPB output y.txt           (one yᵢ per line)
     N_pts        new interior subdivisions per flagged interval (default: 10)
     min_width    minimum interval width to bother refining  (default: 1e-6)
     out_sp_file  where to write the augmented sample list
                  (default: overwrites sp_file in place)

   EXIT CODES
     0  converged — no refinable negative intervals remain, out_sp_file unchanged
     1  refined   — new points written to out_sp_file, loop continues
     2  error     — bad arguments or unreadable files
   ================================================================ *)


(* ----------------------------------------------------------------
   1.  ARGUMENT PARSING
       $ScriptCommandLine = {scriptname, arg1, arg2, ...}
       Works identically for  wolframscript -file  and  math -script.
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
minWidth  = If[Length[myArgs] >= 4, ToExpression[myArgs[[4]]], 10^-6];
outFile   = If[Length[myArgs] >= 5, myArgs[[5]], spFile];

Print["=== negative_region.m ==="];
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

(* Fix C/Fortran-style scientific notation: e.g. "1.5e+02" → "1.5*^+02"
   Mathematica does not recognise lowercase 'e' or uppercase 'E' as an
   exponent marker; it treats them as the symbol e or Euler's E.
   This regex converts  e±digits  /  E±digits  to  *^±digits  before
   ToExpression is called. *)
fixSciNotation[s_String] := StringReplace[s,
  RegularExpression["[eE]([+-]?\\d+)"] :> "*^$1"
];

(* Parse and sort; keep high precision *)
samplePoints = Sort[SetPrecision[ToExpression /@ (fixSciNotation /@ spRaw), 50]];

Print["Loaded ", Length[samplePoints], " sampling points:"];
Print["  ", samplePoints];
Print[""];


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

yVec = SetPrecision[ToExpression /@ (fixSciNotation /@ yRaw), 50];

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


(* ================================================================
   PROBLEM-SPECIFIC SECTION  ← edit here for your problem
   ----------------------------------------------------------------
   Define fList = {f1, f2, f3} to match g3_ExtremalEFT_2.m.
   Each function accepts TWO arguments: a numeric x and an integer J.
   No symbolic form is required — only numeric evaluation is used.

   Jmax, Jlist:  the exact discrete spin values to check at each
                 x-midpoint.  Must match g3_ExtremalEFT_2.m.

   xLeft, xRight: the physical domain boundaries.  The check covers
                  [xLeft, xRight] in full, including the boundary
                  intervals [xLeft, x_min] and [x_max, xRight] that
                  lie outside the current set of sample points.
                  SDPB only enforces positivity AT the sample points;
                  these outer regions are invisible to it otherwise.
   ================================================================ *)
maVal = SetPrecision[0.150, 50];

(* dispersion representation of Wilson coefficients *)
(* All functions now precompute sp and mA numerically with N[...,50] *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[1/2*Sqrt[sp/(sp - 4*mA^2)] * (sp^(-3) + (-4*mA^2 + sp)^(-3)), 50]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(-Sqrt[sp/(sp - 4*mA^2)] * ((-3 + J*(1 + J)*(-4*mA^2 + sp)^3*(sp^(-3) + (-4*mA^2 + sp)^(-3)))/(-4*mA^2 + sp)^4)), 50]
];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(4*mA - sp)^(-5) + (4 - (-2 + J)*J*(1 + J)*(3 + J))/(4*sp^5) + (2*J*(1 + J))/(sp*(-4*mA + sp)^4) - ((-1 + J)*J*(1 + J)*(2 + J))/(4*sp^2*(-4*mA + sp)^3), 50]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 50]
];

X62[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(-2304*mA^4 + 1152*mA^2*sp + (-72 + J*(1 + J)*(-18 + J + J^2))*sp^2))/sp^5))/(576*(-4*mA^2 + sp)^7), 50]
];

X72[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(-47 + J + J^2))/(-4*mA^2 + sp)^8) + ((-1 + J)*(2 + J)*(((-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*sp^3)/(4*mA^2 - sp)^5 + 3600/(-4*mA^2 + sp)^2))/sp^6))/14400, 50]
];

X82[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(((-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(-38 + J + J^2))/(4*mA^2 - sp)^9 + ((-1 + J)*(2 + J)*(-(((-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*sp^4)/(-4*mA^2 + sp)^6) + 129600/(-4*mA^2 + sp)^2))/sp^7))/518400, 50]
];

X92[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(2754 + J*(1 + J)*(-119 + J + J^2)))/(-4*mA^2 + sp)^10) + ((-1 + J)*(2 + J)*(((-6 + J)*(-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*(7 + J)*sp^5)/(4*mA^2 - sp)^7 + 6350400/(-4*mA^2 + sp)^2))/sp^8))/25401600, 50]
];

X102[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-10 + J)*(-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(11 + J)*(2232 + (-10 + J)*J*(1 + J)*(11 + J)))/(-4*mA^2 + sp)^11) + ((-1 + J)*(2 + J)*(-(((-7 + J)*(-6 + J)*(-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*(7 + J)*(8 + J)*sp^6)/(-4*mA^2 + sp)^8) + 406425600/(-4*mA^2 + sp)^2))/sp^9))/1625702400, 50]
];

(* Large J limit *)
LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-(-4*mA^2 + sp)^(-11) - 1/(sp^3*(-4*mA^2 + sp)^8)))/1625702400, 50]
];

Jmax = 40;
Jlist = Range[0, Jmax, 2];

fList = {g20, g31, n4, X52, X62, X72, X82, X92, X102};

(* large J limit *)
(* 0& is a constant function of 0 *)

extraTriplet = {0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, LargeJ};

xLeft  = 0;   (* physical domain left endpoint  — check includes [xLeft,  x_min] *)
xRight = 1;   (* physical domain right endpoint — check includes [x_max, xRight] *)

(* ================================================================
   END OF PROBLEM-SPECIFIC SECTION
   ================================================================ *)


(* --- Dimensional consistency check ---
   y.txt must have exactly as many lines as there are functions in
   fList.  A mismatch means the wrong y.txt or fList was supplied. *)
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



(* F(x, J) = sum_{k=1}^{n} fList[[k]][x, J] * yVec[[k]]
   Evaluated for each (midpoint, J) pair in the interval check below. *)
F[x_?NumericQ, J_?IntegerQ] := Sum[yVec[[k]] * fList[[k]][x, J], {k, Length[yVec]}];

(* large J limit check *)
X[x_?NumericQ] := Sum[yVec[[k]] * extraTriplet[[k]][x], {k, Length[yVec]}];

(* Safety: clamp midpoints near known singularity at x -> 1 (sp = 1/(1-x))
   and provide a robust numeric evaluator that returns $Failed on
   non-finite or exceptional results. *)
singularTol = 10^-12;   (* distance from x=1 to avoid sp singularity *)
safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  val = Quiet[Check[N[F[x0, J]], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];

(* Safe evaluator for large-J functional X[x] *)
safeX[x_?NumericQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  val = Quiet[Check[N[X[x0]], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];

(* ----------------------------------------------------------------
   5.  MIDPOINT CHECK OVER CONSECUTIVE INTERVALS
       For each pair (samplePoints[[i]], samplePoints[[i+1]]):
         - compute midpoint mid and interval width w = xb - xa
         - evaluate F(mid, J) for EVERY J in Jlist
         - fm = Min over J  (worst-case constraint violation)
         - if fm < 0 AND w >= minWidth: flag for x-refinement
         - if fm < 0 AND w <  minWidth: negative but below threshold
   ---------------------------------------------------------------- *)


nIntervals       = Length[samplePoints] - 1;
(* ----------------------------------------------------------------
   Build the complete interval list:
     LEFT boundary  [xLeft,            samplePoints[[1]]  ] if x_min > xLeft
     INTERIOR pairs [samplePoints[[i]], samplePoints[[i+1]]] for i = 1..n-1
     RIGHT boundary [samplePoints[[-1]], xRight            ] if x_max < xRight

   SDPB guarantees F(xi, J) >= 0 at every sample point xi, so the
   interior pairs are automatically safe AT their endpoints.  The
   boundary intervals [xLeft, x_min] and [x_max, xRight] are NEVER
   covered by any sample point — they are invisible to SDPB and must
   be probed explicitly here.  Omitting them was the root cause of
   missing the negative region at x ∈ [0, x_min] seen in z.txt.
   ---------------------------------------------------------------- *)
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
flaggedIntervals = {};   (* intervals to refine: min_J F(mid,J)<0 AND width>=minWidth *)
tinyNegative     = {};   (* intervals with min_J F(mid,J)<0 but width<minWidth        *)
midpointValues   = {};   (* {mid, min_J F(mid,J)} for each interval, for diagnostics  *)

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

  (* Evaluate F(mid, J) safely for every J in Jlist; handle failures gracefully *)
  fmByJ = Table[safeF[mid, Jlist[[j]]], {j, Length[Jlist]}];
  failedCount = Count[fmByJ, $Failed];
  If[failedCount > 0,
    Print["  WARNING: ", failedCount,
          " non-finite F(mid,J) evaluation(s) at mid=", mid,
          " — proceeding with available J values (if any)."]
  ];

  (* evaluate large-J functional X at midpoint *)
  xVal = safeX[mid];
  If[xVal === $Failed,
    Print["  NOTE: X(mid) evaluation failed at mid=", mid]
  ];

  validPairs = Select[Transpose[{Jlist, fmByJ}], Last[#] =!= $Failed &];

  If[validPairs == {} && xVal === $Failed,
    Print["  ERROR: F(mid,J) failed for all J and X(mid) failed at mid=", mid, "; flagging interval for refinement."];
    AppendTo[flaggedIntervals, {xa, xb}];
    AppendTo[midpointValues, {mid, $Failed, "all_failed"}];
    Continue[];
  ];

  (* compute minima from available sources *)
  If[validPairs == {},
    fmJ = $Failed,
    vals = Last /@ validPairs;
    js = First /@ validPairs;
    fmJ = Min[vals];
    worstJIndex = First[Ordering[vals, 1]];
    worstJcandidate = js[[worstJIndex]];
  ];

  (* combine fmJ and xVal to get overall worst-case value and source *)
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
            "  min value=", fm, "  (source=", worstSource, If[worstSource=="J", ", worst J=", ""],
            If[worstSource=="J", worstJ, ""], ")",
            "  *** NEGATIVE, will refine ***"]
  ],
  {i, nAllIntervals}
];
Print[""];


(* ----------------------------------------------------------------
   6.  CONVERGENCE CHECK
       Converged when there are no refinable intervals, i.e. every
       interval with F(mid)<0 is already narrower than minWidth.
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
       For each flagged [xa, xb] (from observation.md):
         x*  = (xa + xb) / 2
         s   = (xb - xa) / nPts
         new = { x* + (k - nPts/2) * s  :  k = 0, ..., nPts }
             = { xa, xa+s, ..., x*, ..., xb-s, xb }
       This places nPts+1 points spanning [xa, xb] exactly, including
       the endpoints xa and xb.  In section 8 these will be merged
       with the original samplePoints; the dedup step removes any
       exact coincidences between the two sets.
   ---------------------------------------------------------------- *)

newPoints = Flatten[
  Table[
    Module[{xa, xb, xStar, s},
      xa    = interval[[1]];
      xb    = interval[[2]];
      (* Clamp boundary intervals away from singularities at x=0 and x=1 *)
      If[Abs[xa - xLeft] < singularTol, xa = xLeft + singularTol];
      If[Abs[xb - xRight] < singularTol, xb = xRight - singularTol];
      xStar = (xa + xb) / 2;
      s     = (xb - xa) / nPts;
      Print["  Flagged interval [", interval[[1]], ", ", interval[[2]], "] (clamped to [", xa, ", ", xb, "]):"];
      Print["    x* = ", xStar, "  s = ", s];
      Table[SetPrecision[xStar + (k - nPts/2) * s, 50], {k, 0, nPts}]
    ],
    {interval, flaggedIntervals}
  ]
];

Print[""];
Print["New candidate points (before dedup): ", Length[newPoints]];

(* deduplication tolerance (absolute) — configurable *)
dedupTol = 10^-10;


(* ----------------------------------------------------------------
   8.  MERGE ORIGINAL AND NEW POINTS, THEN DEDUPLICATE AND SORT
       ----------------------------------------------------------------
   WHY WE KEEP THE ORIGINAL SAMPLING POINTS (observation.md)

   Discarding the original samplePoints and writing only newPoints
   causes oscillatory instability across iterations:

     Iter 1: original points cover [x_min, x_max].  Negative region
             found in [xLeft, x_min].  newPoints are all near x_min.
     Iter 2 (old code): sampling_points.txt = newPoints only.
             SDPB is constrained only near x_min; the interval
             [x'_max, xRight] is now completely unconstrained.
             SDPB exploits this gap → negative region appears on
             the right side.
     Iter 3: newPoints now cover the right side.  Left side is
             unconstrained again → negative region jumps back left.
     → The negative region oscillates between the two boundaries
       and never converges.  (Observed in practice, noted in
       observation.md under "Unexpected Instability".)

   FIX: always include the original samplePoints in the output.
   They act as a GLOBAL SKELETON that prevents large unconstrained
   gaps from opening on either side of the locally refined region.
   The adaptive new points provide LOCAL refinement near each
   negative region.  Together they ensure both coverage and
   convergence.  The extra cost (more SDPB blocks per iteration)
   is acceptable and necessary for stability.
   ----------------------------------------------------------------
   Dedup tolerance 10^-10 removes exact duplicates (e.g. an
   endpoint of a flagged interval that coincides with an existing
   sample point to full precision).  Distinct points that happen
   to be close but numerically different are preserved.
   ---------------------------------------------------------------- *)

allPoints = Sort[Join[samplePoints, newPoints]];

dedupPoints = Fold[
  Function[{kept, pt},
    If[Abs[pt - Last[kept]] < dedupTol, kept, Append[kept, pt]]
  ],
  {First[allPoints]},
  Rest[allPoints]
];

(* Safety filter: remove any points exactly at x=0 or x=1 to avoid singularities *)
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
       Original samplePoints + new refined points, deduplicated.
       One high-precision number per line, plain decimal notation
       (no scientific notation) with 50 significant digits.
   ---------------------------------------------------------------- *)

(* Format a number as a plain decimal string — NEVER scientific notation.
   Uses RealDigits to decompose the number, then manually assembles the
   string with leading zeros or integer part as needed.
     0.001      → "0.0010000000000000000000000000000000000000000000000000"
     1.0*10^-12 → "0.0000000000010000000000000000000000000000000000000000"
   This avoids ToString[N[...]] which renders exponents as multi-line
   superscripts in plain text output. *)
formatPlainDecimal[x_?NumericQ] := Module[
  {ax, sign, digits, dpos, ndig = 50, dstr, result},
  If[x == 0, Return["0." <> StringJoin[Table["0", {ndig}]]]];
  sign = If[Negative[x], "-", ""];
  ax   = SetPrecision[Abs[x], ndig];
  {digits, dpos} = RealDigits[ax, 10, ndig];
  (* RealDigits can return negative entries at the tail for
     numbers that don't fill ndig digits; replace with 0 *)
  dstr = StringJoin[ToString /@ Replace[digits, n_ /; n < 0 :> 0, {1}]];
  result = Which[
    dpos >= ndig,
      (* All digits are before the decimal point *)
      dstr <> StringJoin[Table["0", {dpos - ndig}]] <> ".0",
    dpos > 0,
      (* Decimal point falls within the digit string *)
      StringTake[dstr, dpos] <> "." <> StringDrop[dstr, dpos],
    dpos == 0,
      (* 0.d1d2d3... *)
      "0." <> dstr,
    dpos < 0,
      (* 0.000...d1d2d3... — need -dpos leading zeros after "0." *)
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

Quit[1];   (* exit code 1 = adaptive loop should continue *)
