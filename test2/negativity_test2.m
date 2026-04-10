(* ================================================================
   negative_region.m
   ----------------------------------------------------------------
   PURPOSE
     After each SDPB run, read the current sampling points and the
     solver's y.txt, reconstruct the functional

         F(x, J) = f₁(x,J)·y₁ + f₂(x,J)·y₂ + f₃(x,J)·y₃

     for EVERY J in Jlist (discrete spin, constrained exactly), and
     identify x-intervals where F is negative at the midpoint for
     any J.  Also checks the x-independent J→∞ limit constraint.

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
     The extraTriplet = {0,0,2} encodes the J→∞ limit constraint
     0·y₁ + 0·y₂ + 2·y₃ ≥ 0.  This is x-independent and is
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

(* Parse and sort; keep high precision *)
samplePoints = Sort[SetPrecision[ToExpression /@ spRaw, 50]];

Print["Loaded ", Length[samplePoints], " sampling points:"];
Print["  ", samplePoints];
Print[""];


(* ----------------------------------------------------------------
   3.  READ y.txt
       SDPB writes all n components of y, one per line; blank lines
       and lines starting with "#" are skipped.
       The length of yVec must equal the number of functions in fVec.
       A mismatch is caught by the dimensional check below.
   ---------------------------------------------------------------- *)

If[!FileExistsQ[yFile],
  Print["ERROR: y.txt not found: ", yFile]; Quit[2]];

yRaw = Select[
  ReadList[yFile, String],
  StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &
];

If[Length[yRaw] == 0,
  Print["ERROR: y.txt is empty: ", yFile]; Quit[2]];

yVec = SetPrecision[ToExpression /@ yRaw, 50];

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

   extraTriplet: the J→∞ limiting coefficient vector {c1, c2, c3}
                 such that sum_k ck * yk >= 0 is the large-J
                 constraint.  x-independent; checked once separately.

   xLeft, xRight: the physical domain boundaries.  The check covers
                  [xLeft, xRight] in full, including the boundary
                  intervals [xLeft, x_min] and [x_max, xRight] that
                  lie outside the current set of sample points.
                  SDPB only enforces positivity AT the sample points;
                  these outer regions are invisible to it otherwise.
   ================================================================ *)

f1[x_?NumericQ, J_?IntegerQ] := (1 + x)^2;
f2[x_?NumericQ, J_?IntegerQ] := (1 + x) * (3 - 2*J*(J + 1));
f3[x_?NumericQ, J_?IntegerQ] := 2 * J * (J + 1) * (J*(J + 1) - 8);

fList        = {f1, f2, f3};        (* ← must match fList in g3_ExtremalEFT_2.m *)
Jmax         = 40;
Jlist        = Range[0, Jmax, 2];   (* J = 0, 2, 4, …, 40 — exact constraints  *)
extraTriplet = {0, 0, 2};           (* J→∞ limit: {c1,c2,c3} with ck = lim fk/J^4 *)

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


(* ----------------------------------------------------------------
   5.  MIDPOINT CHECK OVER CONSECUTIVE INTERVALS
       For each pair (samplePoints[[i]], samplePoints[[i+1]]):
         - compute midpoint mid and interval width w = xb - xa
         - evaluate F(mid, J) for EVERY J in Jlist
         - fm = Min over J  (worst-case constraint violation)
         - if fm < 0 AND w >= minWidth: flag for x-refinement
         - if fm < 0 AND w <  minWidth: negative but below threshold

       Before the x-loop, check the x-independent extraTriplet
       constraint once.  A violation there cannot be fixed by
       adding more x sample points.
   ---------------------------------------------------------------- *)

(* --- Check x-independent J→∞ constraint first --- *)
extraCheck = extraTriplet . yVec;   (* dot product: sum_k extraTriplet[[k]]*yVec[[k]] *)
Print["Checking J\[Rule]\[Infinity] (extraTriplet) constraint:"];
If[extraCheck < 0,
  Print["  extraTriplet \[CenterDot] y = ", extraCheck,
        "  *** VIOLATED — x-refinement cannot fix this ***"],
  Print["  extraTriplet \[CenterDot] y = ", extraCheck, "  ok"]
];
Print[""];

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

  (* Evaluate F(mid, J) for every J in Jlist; record the worst (minimum) value *)
  fmByJ  = Table[F[mid, Jlist[[j]]], {j, Length[Jlist]}];
  fm     = Min[fmByJ];
  worstJ = Jlist[[ First[Ordering[fmByJ, 1]] ]];   (* J that achieves the minimum *)
  AppendTo[midpointValues, {mid, fm}];

  Which[
    fm >= 0,
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  min_J F(mid,J)=", fm, "  ok"],

    fm < 0 && width < minWidth,
      AppendTo[tinyNegative, {xa, xb}];
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  min_J F(mid,J)=", fm, "  (worst J=", worstJ, ")",
            "  negative but w < ", minWidth, " — below threshold, skipping"],

    fm < 0 && width >= minWidth,
      AppendTo[flaggedIntervals, {xa, xb}];
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  min_J F(mid,J)=", fm, "  (worst J=", worstJ, ")",
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
       the endpoints — which are kept since the old file is overwritten.
   ---------------------------------------------------------------- *)

newPoints = Flatten[
  Table[
    Module[{xa, xb, xStar, s},
      xa    = interval[[1]];
      xb    = interval[[2]];
      xStar = (xa + xb) / 2;
      s     = (xb - xa) / nPts;
      Print["  Flagged interval [", xa, ", ", xb, "]:"];
      Print["    x* = ", xStar, "  s = ", s];
      Table[SetPrecision[xStar + (k - nPts/2) * s, 50], {k, 0, nPts}]
    ],
    {interval, flaggedIntervals}
  ]
];

Print[""];
Print["New candidate points (before dedup): ", Length[newPoints]];


(* ----------------------------------------------------------------
   8.  DEDUPLICATE AND SORT NEW POINTS ONLY
       The new points already span every flagged interval, including
       its endpoints (xa, xb). Those endpoints came from the old
       sampling_points.txt, but the old file is being OVERWRITTEN, so
       we only need the new interior grid — not the previous points.
       Dedup tolerance 10^-10 guards against near-identical endpoints
       produced by two adjacent flagged intervals sharing a boundary.
   ---------------------------------------------------------------- *)

allPoints = Sort[newPoints];

dedupPoints = Fold[
  Function[{kept, pt},
    If[Abs[pt - Last[kept]] < 10^-10, kept, Append[kept, pt]]
  ],
  {First[allPoints]},
  Rest[allPoints]
];

Print["New points before dedup : ", Length[allPoints]];
Print["New points after  dedup : ", Length[dedupPoints]];
Print[""];


(* ----------------------------------------------------------------
   9.  OVERWRITE sampling_points.txt WITH THE NEW POINTS ONLY
       One high-precision number per line, 50 significant digits.
   ---------------------------------------------------------------- *)

outLines = StringRiffle[
  ToString[NumberForm[#, {50, 45}]] & /@ dedupPoints,
  "\n"
];

Export[outFile, outLines, "Text"];

Print["Overwrote ", outFile, " with ", Length[dedupPoints], " new sampling points."];
Print[""];
Print["New sampling points:"];
Print["  ", dedupPoints];
Print[""];
Print["STATUS: Refinement needed. Re-run pmp2sdp + sdpb with updated ", outFile, "."];

Quit[1];   (* exit code 1 = adaptive loop should continue *)