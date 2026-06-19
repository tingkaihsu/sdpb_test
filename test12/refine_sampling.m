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
minWidth  = If[Length[myArgs] >= 4, ToExpression[myArgs[[4]]], SetPrecision[1/10000000000, 600]];
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

parseSDPBReal[str_String] := Module[{s = StringTrim[str], expr},
  If[s === "" || StringStartsQ[s, "#"], Return[$Failed]];
(* BUG FIX 1: z.txt (y.txt) writes Euler's number as lowercase 'e'.
   Mathematica's ToExpression treats bare 'e' as a free symbol, not as
   the built-in constant E = 2.71828..., making yVec[[4]] symbolic and
   causing all J != 0 evaluations of F[x,J] and every evaluation of X[x]
   to fail.  Replace every standalone 'e' (not part of an identifier,
   number literal, or backtick precision annotation) with 'E' first. *)
(* yRawFixed = StringReplace[#,
  RegularExpression["(?<![A-Za-z0-9`\\$_])e(?![A-Za-z0-9`_])"] -> "*10^"
] & /@ yRaw;
yVec = SetPrecision[ToExpression /@ yRawFixed, 600]; *)

s = StringReplace[
  s,
  RegularExpression["(?<=\\d)[eE]([+-]?\\d+)"] -> "*^$1"
];

expr = Quiet[Check[ToExpression[s], $Failed]];
If[expr === $Failed || !NumericQ[expr],
  $Failed,
  SetPrecision[N[expr, 600], 600]
]
];

yVec = parseSDPBReal /@ yRaw;

(* Guard: if any component is still not purely numeric, abort early with
   a clear message rather than silently producing $Failed in every safeF/safeX call. *)
Do[
  If[!NumericQ[yVec[[k]]],
    Print["ERROR: y.txt component ", k, " is not a number after parsing: ", yVec[[k]]];
    Print["  Raw line was: ", yRaw[[k]]];
    Quit[2]
  ],
  {k, Length[yVec]}
];

Print["y vector (", Length[yVec], " component(s)): ", yVec];
Print[""];


(* ================================================================
   PROBLEM-SPECIFIC SECTION  ← edit here for your problem
   ----------------------------------------------------------------
   Define fVec = {f1, f2, …, fn} to match the functions used in
   Tests2.m.  Each function accepts a numeric argument and returns
   a numeric value.  No symbolic form is required.
   ================================================================ *)
Print["Mass scales are normalized by the first isolated state..."];
Print[""];
J1 = 4;
Print["J1 = ", J1];

m2 = N[1, 600];
Print["m2^2 = ", m2];
J2 = 2;
Print["J2 = ", J2];

mgap = N[3, 600];

maVal = N[1/1000, 600];

Print["m_gap^2  = ", mgap];
Print["mA = ", maVal];

(* forward limit: use our own convention *)
g2shft[x_?NumericQ] := Module[
  {sp, mA},
  sp = SetPrecision[m2, 600];
  mA = SetPrecision[maVal, 600];
  N[(2 Sqrt[sp/(-4 mA^2 + sp)])/(sp - 2 mA^2)^3, 600]
];

n4AAAAshft[x_?NumericQ] := Module[
  {sp, mA, J},
  sp = SetPrecision[m2, 600];
  mA = SetPrecision[maVal, 600];
  J  = J2;
  Re[N[-((243 Sqrt[-(sp/(4 mA^2-sp))] (-6 I (8 mA^3-3 mA sp) LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+Sqrt[-8 mA^2+3 sp] (-4 mA^2+3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (-8 mA^2+3 sp)^(3/2))), 600] ]
];

X52AAAAshft[x_?NumericQ] := Module[
  {sp, mA, J},
  sp = SetPrecision[m2, 600];
  mA = SetPrecision[maVal, 600];
  J  = J2;
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

X62AAAAshft[x_?NumericQ] := Module[
  {sp, mA, J},
  sp = SetPrecision[m2, 600];
  mA = SetPrecision[maVal, 600];
  J  = J2;
  N[((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7)), 600]
];

X72AAAAshft[x_?NumericQ] := Module[
  {sp, mA, J},
  sp = SetPrecision[m2, 600];
  mA = SetPrecision[maVal, 600];
  J  = J2;
  N[((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400), 600]
];

(* g2 > 0 *)
g2[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[mgap/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[((2 Sqrt[sp/(-4 mA^2+sp)])/(sp - 2 mA^2)^3), 600]
];

n4AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp   = N[mgap/(1 - x), 600];
  mA   = N[maVal, 600];
  result = -((243 Sqrt[-(sp/(4 mA^2-sp))] (-6 I (8 mA^3-3 mA sp) LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+Sqrt[-8 mA^2+3 sp] (-4 mA^2+3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (-8 mA^2+3 sp)^(3/2)));
  Re[N[result, 600] ]
];


X52AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[maVal, 600];
  result = (J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6);
  Re[N[result, 600] ]
];

X62AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[mgap/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7)), 600]
];


X72AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[maVal, 600];
  result = ((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400);
  Re[N[result, 600] ]
];
  
(* k = 7, q = 2 *)
LargeJAAAA[x_?NumericQ] := Module[{sp, mA},
  sp = N[mgap/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[-(((1/(-4 mA^2+sp))^(17/2) (-32 mA^6+24 mA^4 sp-6 mA^2 sp^2+sp^3))/(7200 sp^(5/2))), 600]
];

Jmax = 66;
(* FIX: Jlist must start from J=0 to include all even spins in UV continuum.
   The paper imposes positivity over all even spins l in [0, 500].
   Missing J=0,2,4 removes nontrivial positivity conditions. *)
Jlist = Range[0, Jmax, 2];

(* 2g[2,2] - g[2,1] *)
M0[x_?NumericQ,J_?IntegerQ] := {
	{g2[x, J],0,0},
	{0,0,0},
	{0,0,0}
};

M0shft[x_?NumericQ] := {
  {g2shft[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M41[x_?NumericQ, J_?IntegerQ] := {
  {n4AAAA[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M41shft[x_?NumericQ] := {
  {n4AAAAshft[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M51[x_?NumericQ, J_?IntegerQ] := {
  {X52AAAA[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M51shft[x_?NumericQ] := {
  {X52AAAAshft[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M61[x_?NumericQ, J_?IntegerQ] := {
  {X62AAAA[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M61shft[x_?NumericQ] := {
  {X62AAAAshft[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M71[x_?NumericQ, J_?IntegerQ] := {
  {X72AAAA[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M71shft[x_?NumericQ] := {
  {X72AAAAshft[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

M7j1[x_?NumericQ] :={
	{LargeJAAAA[x],0,0},
	{0,0,0},
	{0,0,0}
};

(* FIX: Define LargeJAAAAshft — was referenced in M7j1shft but never defined *)
LargeJAAAAshft[x_?NumericQ] := Module[{sp, mA},
  sp = SetPrecision[m2, 600];
  mA = SetPrecision[maVal, 600];
  N[-(((1/(-4 mA^2+sp))^(17/2) (-32 mA^6+24 mA^4 sp-6 mA^2 sp^2+sp^3))/(7200 sp^(5/2))), 600]
];

M7j1shft[x_?NumericQ] := {
  {LargeJAAAAshft[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

(* null *)
N0[x_?NumericQ] :={
	{0,0,0},
	{0,0,0},
	{0,0,0}
};

f11List ={
	Function[{x,J}, M0[x,J][[1,1]]],
	Function[{x,J}, M41[x,J][[1,1]]],
  Function[{x,J}, M51[x,J][[1,1]]],
  Function[{x,J}, M61[x,J][[1,1]]],
  Function[{x,J}, M71[x,J][[1,1]]]
};

f22List ={
	Function[{x,J}, M0[x,J][[2,2]]],
	Function[{x,J}, M41[x,J][[2,2]]],
  Function[{x,J}, M51[x,J][[2,2]]],
  Function[{x,J}, M61[x,J][[2,2]]],
  Function[{x,J}, M71[x,J][[2,2]]]
};

f33List = {
	Function[{x,J}, M0[x,J][[3,3]]],
	Function[{x,J}, M41[x,J][[3,3]]],
  Function[{x,J}, M51[x,J][[3,3]]],
  Function[{x,J}, M61[x,J][[3,3]]],
  Function[{x,J}, M71[x,J][[3,3]]]
};

f12List ={
	Function[{x,J}, M0[x,J][[1,2]]],
	Function[{x,J}, M41[x,J][[1,2]]],
  Function[{x,J}, M51[x,J][[1,2]]],
  Function[{x,J}, M61[x,J][[1,2]]],
  Function[{x,J}, M71[x,J][[1,2]]]
};

f21List = f12List;

f13List = {
	Function[{x,J}, M0[x,J][[1,3]]],
	Function[{x,J}, M41[x,J][[1,3]]],
  Function[{x,J}, M51[x,J][[1,3]]],
  Function[{x,J}, M61[x,J][[1,3]]],
  Function[{x,J}, M71[x,J][[1,3]]]
};
f31List = f13List;

f23List = {
	Function[{x,J}, M0[x,J][[2,3]]],
	Function[{x,J}, M41[x,J][[2,3]]],
  Function[{x,J}, M51[x,J][[2,3]]],
  Function[{x,J}, M61[x,J][[2,3]]],
  Function[{x,J}, M71[x,J][[2,3]]]
};

f32List = f23List;

f11ShftList = {
  Function[{x}, M0shft[x][[1, 1]]],
  Function[{x}, M41shft[x][[1, 1]]],
  Function[{x}, M51shft[x][[1, 1]]],
  Function[{x}, M61shft[x][[1, 1]]],
  Function[{x}, M71shft[x][[1, 1]]]
};

f22ShftList = {
  Function[{x}, M0shft[x][[2, 2]]],
  Function[{x}, M41shft[x][[2, 2]]],
  Function[{x}, M51shft[x][[2, 2]]],
  Function[{x}, M61shft[x][[2, 2]]],
  Function[{x}, M71shft[x][[2, 2]]]
};

f33ShftList = {
  Function[{x}, M0shft[x][[3, 3]]],
  Function[{x}, M41shft[x][[3, 3]]],
  Function[{x}, M51shft[x][[3, 3]]],
  Function[{x}, M61shft[x][[3, 3]]],
  Function[{x}, M71shft[x][[3, 3]]]
};

f12ShftList = {
  Function[{x}, M0shft[x][[1, 2]]],
  Function[{x}, M41shft[x][[1, 2]]],
  Function[{x}, M51shft[x][[1, 2]]],
  Function[{x}, M61shft[x][[1, 2]]],
  Function[{x}, M71shft[x][[1, 2]]]
};

f21ShftList = f12ShftList;

f13ShftList = {
  Function[{x}, M0shft[x][[1, 3]]],
  Function[{x}, M41shft[x][[1, 3]]],
  Function[{x}, M51shft[x][[1, 3]]],
  Function[{x}, M61shft[x][[1, 3]]],
  Function[{x}, M71shft[x][[1, 3]]]
};

f31ShftList = f13ShftList;

f23ShftList = {
  Function[{x}, M0shft[x][[2, 3]]],
  Function[{x}, M41shft[x][[2, 3]]],
  Function[{x}, M51shft[x][[2, 3]]],
  Function[{x}, M61shft[x][[2, 3]]],
  Function[{x}, M71shft[x][[2, 3]]]
};

f32ShftList = f23ShftList;

j11List = {
	Function[{x}, N0[x][[1,1]]],
	Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, M7j1[x][[1,1]]]
};

j22List = {
	Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, M7j1[x][[2,2]]]
};

j33List = {
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, M7j1[x][[3,3]]]
};

j12List = {
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, M7j1[x][[1,2]]]
};
j21List = j12List;

j13List = {
  Function[{x}, N0[x][[1,3]]],
	Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
  Function[{x}, M7j1[x][[1,3]]]
};
j31List = j13List;

j23List = {
  Function[{x}, N0[x][[2,3]]],
	Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
  Function[{x}, M7j1[x][[2,3]]]
};
j32List = j23List;

G2[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[((2 Sqrt[sp/(-4 mA^2+sp)])/(sp - 2 mA^2)^3), 600]
];

N4AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp   = N[1/(1 - x), 600];
  mA   = N[maVal, 600];
  result = -((243 Sqrt[-(sp/(4 mA^2-sp))] (-6 I (8 mA^3-3 mA sp) LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+Sqrt[-8 mA^2+3 sp] (-4 mA^2+3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (-8 mA^2+3 sp)^(3/2)));
  Re[N[result, 600] ]
];

x52AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = (J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6);
  Re[N[result, 600] ]
];

x62AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  (* FIX: was copy-pasted from x52AAAA; now uses correct X62 formula *)
  N[((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7)), 600]
];

x72AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[(J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400, 600]
];

m1 = SetPrecision[1/3, 600];

norm = {G2[m1, J1], N4AAAA[m1, J1], x52AAAA[m1, J1], x62AAAA[m1, J1], x72AAAA[m1, J1]};

obj  = {-1, 0, 0, 0, 0};

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
F31[x_?NumericQ, J_?IntegerQ] := F13[x, J];

F23[x_?NumericQ, J_?IntegerQ] := Sum[ yVec[[k]] * f23List[[k]][x, J], {k, Length[yVec]} ];
F32[x_?NumericQ, J_?IntegerQ] := F23[x, J];

J11[x_?NumericQ] := Sum[ yVec[[k]] * j11List[[k]][x], {k, Length[yVec]} ];
J22[x_?NumericQ] := Sum[ yVec[[k]] * j22List[[k]][x], {k, Length[yVec]} ];
J33[x_?NumericQ] := Sum[ yVec[[k]] * j33List[[k]][x], {k, Length[yVec]} ];

J12[x_?NumericQ] := Sum[ yVec[[k]] * j12List[[k]][x], {k, Length[yVec]} ];
J21[x_?NumericQ] := J12[x];

J13[x_?NumericQ] := Sum[ yVec[[k]] * j13List[[k]][x], {k, Length[yVec]} ];
J31[x_?NumericQ] := J13[x];

J23[x_?NumericQ] := Sum[ yVec[[k]] * j23List[[k]][x], {k, Length[yVec]} ];
J32[x_?NumericQ] := J23[x];

F[x_?NumericQ, J_?IntegerQ] := {
    {F11[x, J], F12[x, J], F13[x, J]},
    {F21[x, J], F22[x, J], F23[x, J]},
    {F31[x, J], F32[x, J], F33[x, J]}
};

X[x_?NumericQ] := {
    {J11[x], J12[x], J13[x]},
    {J21[x], J22[x], J23[x]},
    {J31[x], J32[x], J33[x]}
};

(* Safety: clamp midpoints near known singularity at x -> 1 (sp = 1/(1-x))
   and provide a robust numeric evaluator that returns $Failed on
   non-finite or exceptional results. *)
singularTol = SetPrecision[1/10000000000, 600];   (* distance from x=1 to avoid sp singularity *)

(* BUG FIX 2-5 for safeF:
   (2) 'mat' and 'egnVals' are now LOCAL Module variables — they were
       previously global, causing stale values from one call to corrupt the next.
   (3) Both N[F[…]] and Eigenvalues[…] are wrapped in Quiet[Check[…,$Failed]]
       so that SILENT failures (returning unevaluated without a message) are
       caught before reaching Min, not just message-generating ones.
   (4) MatrixQ[mat, NumericQ] pre-flight guard rejects symbolic or partial
       matrices before Eigenvalues is called, giving a clean $Failed path.
   (5) Re /@ egnVals strips tiny numerical imaginary parts that arise from
       floating-point arithmetic on nearly-symmetric matrices; the final
       Im[minEgn] == 0 check ensures no complex value leaks through
       (NumberQ[a+bI] is True, so the old MemberQ guard was insufficient). *)
safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, mat, egnVals, minEgn},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];

  (* Step 1: evaluate the matrix, catching both message and silent failures *)
  mat = Quiet[Check[N[F[x0, J], 600], $Failed]];
  If[mat === $Failed, Return[$Failed]];

  (* Step 2: reject symbolic / non-numeric matrices immediately *)
  If[!MatrixQ[mat, NumericQ], Return[$Failed]];

  (* Step 3: compute eigenvalues, catching failures *)
  egnVals = Quiet[Check[Eigenvalues[mat], $Failed]];
  If[egnVals === $Failed || !VectorQ[egnVals, NumericQ], Return[$Failed]];

  (* Step 4: strip floating-point imaginary noise, then take the minimum *)
  minEgn = Min[Re /@ egnVals];

  (* Step 5: final sanity — result must be a real finite number *)
  If[!NumberQ[minEgn] || Im[minEgn] =!= 0 ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate}, minEgn],
    $Failed,
    minEgn
  ]
];

(* BUG FIX 2-5 for safeX: same set of fixes as safeF above. *)
safeX[x_?NumericQ] := Module[{x0 = x, mat, egnVals, minEgn},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];

  (* Step 1: evaluate the matrix, catching both message and silent failures *)
  mat = Quiet[Check[N[X[x0], 600], $Failed]];
  If[mat === $Failed, Return[$Failed]];

  (* Step 2: reject symbolic / non-numeric matrices immediately *)
  If[!MatrixQ[mat, NumericQ], Return[$Failed]];

  (* Step 3: compute eigenvalues, catching failures *)
  egnVals = Quiet[Check[Eigenvalues[mat], $Failed]];
  If[egnVals === $Failed || !VectorQ[egnVals, NumericQ], Return[$Failed]];

  (* Step 4: strip floating-point imaginary noise, then take the minimum *)
  minEgn = Min[Re /@ egnVals];

  (* Step 5: final sanity — result must be a real finite number *)
  If[!NumberQ[minEgn] || Im[minEgn] =!= 0 ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate}, minEgn],
    $Failed,
    minEgn
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

nIntervals = Length[samplePoints] - 1;

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
      If[samplePoints[[1]] > xLeft  + singularTol, "(active)", "(skipped — x_min = xLeft)"]];
Print["  Right boundary : [", samplePoints[[-1]], ", ", xRight, "]  ",
      If[samplePoints[[-1]] < xRight - singularTol, "(active)", "(skipped — x_max = xRight)"]];
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

  (* Evaluate large-J functional X at midpoint *)
  xVal = safeX[mid];
  If[xVal === $Failed,
    Print["  NOTE: X(mid) evaluation failed at mid=", mid]
  ];

  (* If every source failed, flag conservatively and continue *)
  If[failedCount == Length[Jlist] && xVal === $Failed,
    Print["  Interval [", xa, ", ", xb, "]  w=", width, "  mid=", mid];
    Print["    ERROR: all evaluations failed — flagging for refinement"];
    AppendTo[flaggedIntervals, {xa, xb}];
    AppendTo[midpointValues, {mid, $Failed, 0}];
    Continue[]
  ];

  (* ---------------------------------------------------------------
     Collect every individual violation.
     Each entry: {"J", J_value, min_eigenvalue}  or  {"X", None, min_eigenvalue}
     A violation exists whenever the minimum eigenvalue of F[mid,J]
     or X[mid] is strictly negative — independently of the other sources.
     We do NOT aggregate into a single scalar; every (x,J) and (x,X)
     violation is recorded and reported separately.
     --------------------------------------------------------------- *)
  violations = {};

  (* Check each J value independently *)
  Do[
    With[{val = fmByJ[[jIdx]], Jval = Jlist[[jIdx]]},
    (* < -10^-10 instead of < 0 *)
      If[val =!= $Failed && val < -1*singularTol,
        AppendTo[violations, {"J", Jval, val}]
      ]
    ],
    {jIdx, Length[Jlist]}
  ];

  (* Check the large-J functional independently *)
  If[xVal =!= $Failed && xVal < -1*singularTol,
    AppendTo[violations, {"X", None, xVal}]
  ];

  (* --- Report interval header then every individual violation --- *)
  Print["  Interval [", xa, ", ", xb, "]  w=", width, "  mid=", mid];

  If[Length[violations] > 0,

    (* Report EACH (x, J) or (x, X) violation on its own line *)
    Do[
      With[{viol = violations[[v]]},
        If[viol[[1]] === "J",
          Print["    *** NEGATIVE: x=", mid,
                "  J=", viol[[2]],
                "  min_eigenvalue=", viol[[3]]],
          Print["    *** NEGATIVE: x=", mid,
                "  X (large-J limit)",
                "  min_eigenvalue=", viol[[3]]]
        ]
      ],
      {v, Length[violations]}
    ];

    (* Flag the interval for refinement, or mark as tiny if too narrow *)
    If[width < minWidth,
      AppendTo[tinyNegative, {xa, xb}];
      Print["    -> negative (", Length[violations], " violation(s)) but width=", width,
            " < min_width=", minWidth, " — below threshold, skipping"],
      AppendTo[flaggedIntervals, {xa, xb}];
      Print["    -> will refine: new sampling point at x=", mid,
            "  (", Length[violations], " violation(s))"]
    ],

    (* No violations at this midpoint: all F(mid,J) >= 0 and X(mid) >= 0 *)
    Print["    ok - F(mid,J) >= 0 for all J and X(mid) >= 0"]
  ];

  (* Diagnostic record: worst (most negative) value among violations,
     or smallest positive value when there are none *)
  negVals = Last /@ violations;
  worstVal = If[Length[negVals] > 0,
    Min[negVals],
    posVals = Select[
      Join[DeleteCases[fmByJ, $Failed], If[xVal === $Failed, {}, {xVal}]],
      # >= 0 &];
    If[Length[posVals] > 0, Min[posVals], $Failed]
  ];
  AppendTo[midpointValues, {mid, worstVal, Length[violations]}],

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
dedupTol = SetPrecision[1/10000000000, 600];


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
     1.0*10^-10 → "0.0000000000010000000000000000000000000000000000000000"
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
