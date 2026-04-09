(* ================================================================
   negative_region.m
   ----------------------------------------------------------------
   PURPOSE
     After each SDPB run, read the current sampling points and the
     solver's y.txt, reconstruct the functional

         F(x) = f1(x)*y1 + f2(x)*y2

     and identify intervals between consecutive sampling points where
     F is negative, by evaluating F at the midpoint of each interval.

   WHY MIDPOINTS ONLY (not a fine grid scan)
     SDPB already guarantees F(xᵢ) >= 0 at every sampling point xᵢ.
     Therefore we do NOT need to scan the whole real line. We only
     need to probe each open interval (xᵢ, xᵢ₊₁). A single midpoint
     evaluation per interval is sufficient: if F((xᵢ+xᵢ₊₁)/2) < 0
     the interval is flagged for refinement. This is both faster and
     more principled than a fine-grid sweep.

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
       SDPB writes one component per line; blank lines and lines
       starting with "#" are skipped.
       The normalization n = (1,0) fixes y[[1]] = 1.  y[[2]] is the
       free optimisation variable returned by the solver.
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

(* Guard: if SDPB only wrote the free component, prepend the fixed y1=1 *)
If[Length[yVec] == 1, yVec = Prepend[yVec, 1]];

Print["y vector: ", yVec];
Print[""];


(* ----------------------------------------------------------------
   4.  FUNCTIONAL DEFINITION  ← edit here for your problem
       F(x) = Sum_k  f_k(x) * y_k
   ---------------------------------------------------------------- *)

f1[x_?NumericQ] := 1 + x^4;
f2[x_?NumericQ] := x^4/12 + x^2;

F[x_?NumericQ] := yVec[[1]] * f1[x] + yVec[[2]] * f2[x];


(* ----------------------------------------------------------------
   5.  MIDPOINT CHECK OVER CONSECUTIVE INTERVALS
       For each pair (samplePoints[[i]], samplePoints[[i+1]]):
         - compute midpoint m and interval width w = xb - xa
         - evaluate F(m)
         - if F(m) < 0 AND w >= minWidth: flag for refinement
         - if F(m) < 0 AND w <  minWidth: negative but below the
           stopping threshold — skip, do not generate new points
   ---------------------------------------------------------------- *)

nIntervals       = Length[samplePoints] - 1;
flaggedIntervals = {};   (* intervals to refine: F(mid)<0 AND width>=minWidth *)
tinyNegative     = {};   (* intervals with F(mid)<0 but width<minWidth        *)
midpointValues   = {};   (* for diagnostic printing *)

Print["Checking ", nIntervals, " interval(s) between consecutive sample points:"];
Print["  (stopping threshold: min_width = ", minWidth, ")"];
Print[""];

Do[
  xa    = samplePoints[[i]];
  xb    = samplePoints[[i + 1]];
  width = xb - xa;
  mid   = (xa + xb) / 2;
  fm    = F[mid];
  AppendTo[midpointValues, {mid, fm}];

  Which[
    fm >= 0,
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  F(mid)=", fm, "  ok"],

    fm < 0 && width < minWidth,
      AppendTo[tinyNegative, {xa, xb}];
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  F(mid)=", fm,
            "  negative but w < ", minWidth, " — below threshold, skipping"],

    fm < 0 && width >= minWidth,
      AppendTo[flaggedIntervals, {xa, xb}];
      Print["  [", xa, ", ", xb, "]  w=", width,
            "  F(mid)=", fm, "  *** NEGATIVE, will refine ***"]
  ],
  {i, nIntervals}
];
Print[""];


(* ----------------------------------------------------------------
   6.  CONVERGENCE CHECK
       Converged when there are no refinable intervals, i.e. every
       interval with F(mid)<0 is already narrower than minWidth.
   ---------------------------------------------------------------- *)

If[Length[flaggedIntervals] == 0,
  If[Length[tinyNegative] == 0,
    Print["CONVERGED: F(midpoint) >= 0 for all ", nIntervals, " intervals."],
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