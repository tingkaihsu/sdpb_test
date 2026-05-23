(* ----------------------------------------------------------------
   PURPOSE
     After each SDPB run, read the current sampling points and the
     solver's y.txt, reconstruct the functional

         F(x) = f₁(x)·y₁ + f₂(x)·y₂ + … + fₙ(x)·yₙ

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
samplePoints = Sort[SetPrecision[ToExpression /@ spRaw, 600]];

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

yVec = SetPrecision[ToExpression /@ yRaw, 600];

Print["y vector (", Length[yVec], " component(s)): ", yVec];
Print[""];


(* ================================================================
   PROBLEM-SPECIFIC SECTION  ← edit here for your problem
   ----------------------------------------------------------------
   Define fVec = {f1, f2, …, fn} to match the functions used in
   Tests2.m.  Each function accepts a numeric argument and returns
   a numeric value.  No symbolic form is required.
   ================================================================ *)
maVal = N[1/100, 600];

Print["mA = ", maVal]

(* forward limit: use our own convention *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[-(2*Sqrt[sp/(-4*mA^2 + sp)])/(2*mA^2 - sp)^3, 600]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[-(Sqrt[sp/(-4*mA^2 + sp)]*(-3 - (2*J*(1 + J)*(2*mA^2 - sp))/(-4*mA^2 + sp)))/(-2*mA^2 + sp)^4, 600]
];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp   = N[1/(1 - x), 600];
  mA   = N[maVal, 600];
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
  Re[N[result, 600] ]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

(* Large J limit *)
LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-32*mA^6 + 24*mA^4*sp - 6*mA^2*sp^2 + sp^3))/(18*sp^3*(-4*mA^2 + sp)^6), 600]
];

Jmax = 60;
Jlist = Range[0, Jmax, 2];

(* Matrix *)
M0[x_?NumericQ,J_?IntegerQ] := {
	{g20[x,J],0,0},

	{0,0,0},
	{0,0,0}
};

M1[x_?NumericQ,J_?IntegerQ] := {
	{g31[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M2[x_?NumericQ, J_?IntegerQ] :={
	{n4[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M3[x_?NumericQ, J_?IntegerQ] :={
	{X52[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M4[x_?NumericQ] :={
	{LargeJ[x],0,0},
	{0,0,0},
	{0,0,0}
}

(* null *)
N0[x_?NumericQ] :={
	{0,0,0},
	{0,0,0},
	{0,0,0}
};

f11List ={
	Function[{x,J}, M0[x,J][[1,1]]],
	Function[{x,J}, M1[x,J][[1,1]]],
	Function[{x,J}, M2[x,J][[1,1]]],
	Function[{x,J}, M3[x,J][[1,1]]]
};

f22List ={
	Function[{x,J}, M0[x,J][[2,2]]],
	Function[{x,J}, M1[x,J][[2,2]]],
	Function[{x,J}, M2[x,J][[2,2]]],
	Function[{x,J}, M3[x,J][[2,2]]]
};

f33List = {
	Function[{x,J}, M0[x,J][[3,3]]],
	Function[{x,J}, M1[x,J][[3,3]]],
	Function[{x,J}, M2[x,J][[3,3]]],
	Function[{x,J}, M3[x,J][[3,3]]]
};

f12List ={
	Function[{x,J}, M0[x,J][[1,2]]],
	Function[{x,J}, M1[x,J][[1,2]]],
	Function[{x,J}, M2[x,J][[1,2]]],
	Function[{x,J}, M3[x,J][[1,2]]]
};

f21List = f12List;

f13List = {
	Function[{x,J}, M0[x,J][[1,3]]],
	Function[{x,J}, M1[x,J][[1,3]]],
	Function[{x,J}, M2[x,J][[1,3]]],
	Function[{x,J}, M3[x,J][[1,3]]]
};
f31List = f13List;

f23List = {
	Function[{x,J}, M0[x,J][[2,3]]],
	Function[{x,J}, M1[x,J][[2,3]]],
	Function[{x,J}, M2[x,J][[2,3]]],
	Function[{x,J}, M3[x,J][[2,3]]]
};
f32List = f23List;

j11List = {
	Function[{x}, N0[x][[1,1]]],
	Function[{x}, N0[x][[1,1]]],
	Function[{x}, N0[x][[1,1]]],
	Function[{x}, M4[x][[1,1]]]
};
j22List = {
	Function[{x}, N0[x][[2,2]]],
	Function[{x}, N0[x][[2,2]]],
	Function[{x}, N0[x][[2,2]]],
	Function[{x}, M4[x][[2,2]]]
};
j33List = {
	Function[{x}, N0[x][[3,3]]],
	Function[{x}, N0[x][[3,3]]],
	Function[{x}, N0[x][[3,3]]],
	Function[{x}, M4[x][[3,3]]]
};

j12List = {
	Function[{x}, N0[x][[1,2]]],
	Function[{x}, N0[x][[1,2]]],
	Function[{x}, N0[x][[1,2]]],
	Function[{x}, M4[x][[1,2]]]
};
j21List = j12List;

j13List = {
	Function[{x}, N0[x][[1,3]]],
	Function[{x}, N0[x][[1,3]]],
	Function[{x}, N0[x][[1,3]]],
	Function[{x}, M4[x][[1,3]]]
};
j31List = j13List;

j23List = {
	Function[{x}, N0[x][[2,3]]],
	Function[{x}, N0[x][[2,3]]],
	Function[{x}, N0[x][[2,3]]],
	Function[{x}, M4[x][[2,3]]]
};
j32List = j23List;

(* optimal upper bound *)
norm = {0, -1, 0, 0};
obj  = {-1, 0, 0, 0};

(* ================================================================
   END OF PROBLEM-SPECIFIC SECTION
   ================================================================ *)

xLeft  = SetPrecision[0, 600];   (* physical domain left endpoint  — check includes [xLeft,  x_min] *)
xRight = SetPrecision[1, 600];   (* physical domain right endpoint — check includes [x_max, xRight] *)

(* --- Dimensional consistency check ---
   y.txt must have exactly as many lines as there are functions in
   fVec.  A mismatch means either the wrong y.txt or the wrong fVec
   was supplied; either way the functional would be evaluated
   incorrectly, so we stop immediately. *)
If[Length[yVec] != Length[f11List],
  Print["ERROR: y.txt has ", Length[yVec], " component(s) but fVec has ",
        Length[f11List], " function(s). They must match."];
  Quit[2]
];

(* F(x) = sum_{k=1}^{n} fVec[[k]](x) * yVec[[k]]
   Works for any n >= 1; no change needed when n changes. *)
F11[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f11List[[k]][x, J], {k, Length[yVec]} ];
F22[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f22List[[k]][x, J], {k, Length[yVec]} ];
F33[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f33List[[k]][x, J], {k, Length[yVec]} ];

F12[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f12List[[k]][x, J], {k, Length[yVec]} ];
F21[x_?NumericQ, J_?IntegerQ] := F12[x, J];

F13[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f13List[[k]][x, J], {k, Length[yVec]} ];
F31[x_?NumericQ, J_?IntegerQ] := F31[x, J];

F23[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f23List[[k]][x, J], {k, Length[yVec]} ];
F32[x_?NumericQ, J_?IntegerQ] := F23[x, J];

J11[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f11List[[k]][x, J], {k, Length[yVec]} ];
J22[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f22List[[k]][x, J], {k, Length[yVec]} ];
J33[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f33List[[k]][x, J], {k, Length[yVec]} ];

J12[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f12List[[k]][x, J], {k, Length[yVec]} ];
J21[x_?NumericQ, J_?IntegerQ] := F12[x, J];

J13[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f13List[[k]][x, J], {k, Length[yVec]} ];
J31[x_?NumericQ, J_?IntegerQ] := J31[x, J];

J23[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f23List[[k]][x, J], {k, Length[yVec]} ];
J32[x_?NumericQ, J_?IntegerQ] := J23[x, J];

F[x_?NumericQ, J_?IntegerQ] := {
    {F11[x, J], F12[x, J], F13[x, J]},
    {F21[x, J], F22[x, J], F23[x, J]},
    {F31[x, J], F32[x, J], F33[x, J]}
};

X[x_?NumericQ] := {
    {J11[x, 0], J12[x, 0], J13[x, 0]},
    {J21[x, 0], J22[x, 0], J23[x, 0]},
    {J31[x, 0], J32[x, 0], J33[x, 0]}
};

(* Safety: clamp midpoints near known singularity at x -> 1 (sp = 1/(1-x))
   and provide a robust numeric evaluator that returns $Failed on
   non-finite or exceptional results. *)
singularTol = 10^-12;   (* distance from x=1 to avoid sp singularity *)

safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  temp = N[F[x0, J], 600];
  egnVal = Eigenvalues[temp];
  val = Quiet[Check[Min[egnVal], $Failed] ];

  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];

safeX[x_?NumericQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  temp = N[X[x0], 600];
  egnVal = Eigenvalues[temp];
  val = Quiet[Check[Min[egnVal], $Failed] ];

  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];
(* ----------------------------------------------------------------
   5.  MIDPOINT CHECK OVER CONSECUTIVE INTERVALS
       For each pair (samplePoints[[i]], samplePoints[[i+1]]):
         - compute midpoint m and interval width w = xb - xa
         - evaluate F(m)
         - if F(m) < 0 AND w >= minWidth: flag for refinement
         - if F(m) < 0 AND w <  minWidth: negative but below the
           stopping threshold — skip, do not generate new points
   ---------------------------------------------------------------- *)

interiorPairs = Table[
  {samplePoints[[i]], samplePoints[[i + 1]]},
  {i, nIntervals}
];
allIntervals = Join[
  If[samplePoints[[1]]  > xLeft  + singularTol, {{xLeft,  samplePoints[[1]]}},  {}],
  interiorPairs,
  If[samplePoints[[-1]] < xRight - singularTol, {{samplePoints[[-1]], xRight}}, {}]
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
      Table[SetPrecision[xStar + (k - nPts/2) * s, 600], {k, 0, nPts}]
    ],
    {interval, flaggedIntervals}
  ]
];

Print[""];
Print["New candidate points (before dedup): ", Length[newPoints]];

(* deduplication tolerance (absolute) — configurable *)
dedupTol = 10^-12;


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
       (no scientific notation) with 600 significant digits.
   ---------------------------------------------------------------- *)

(* Format a number as a plain decimal string — NEVER scientific notation.
   Uses RealDigits to decompose the number, then manually assembles the
   string with leading zeros or integer part as needed.
     0.001      → "0.001000000000000000000000000000000000000000000000₀₀₀"
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
