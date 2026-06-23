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
samplePoints = Sort[SetPrecision[ToExpression /@ (fixSciNotation /@ spRaw), 600]];

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

yVec = SetPrecision[ToExpression /@ (fixSciNotation /@ yRaw), 600];

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
(* problem-specific *)
(* 4mA^2 < M^2 *)
m1 = N[1/4, 600];
J1 = 0;
J2 = 2;
mgap = N[166/100, 600];
mAval = N[1/1000, 600];

Print["m1^2 = ", m1];
Print["J1 = ", J1];
Print["J2 = ", J2];
Print["M^2 = ", mgap];
Print["mA = ", mAval];

(* heavy-sum *)
g0[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (2 s)/(-4 mA^4+s^2);
  N[result, 600]
];

x10[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((8 mA^2 + (-2 + J (7 + J)) s)/(2 s^2 (-4 mA^2 + s)));
  N[result, 600]
];

x20[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((-320 mA^4+20 (8+J (7+J)) mA^2 s+(-20+J (7+J) (-13+J (7+J))) s^2)/(20 s^3 (-4 mA^2+s)^2));
  N[result, 600]
];

x30[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((23040 mA^6+1440 (-12+J (7+J)) mA^4 s+36 (120+J (7+J) (-28+J (7+J))) mA^2 s^2+(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) s^3)/(360 s^4 (-4 mA^2+s)^3));
  N[result, 600]
];

x40[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (2580480 mA^8-161280 (16+J (7+J)) mA^6 s-4032 (-240+J (7+J) (-38+J (7+J))) mA^4 s^2-56 (2880+J (7+J) (972+J (7+J) (-62+J (7+J)))) mA^2 s^3+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) s^4)/(10080 s^5 (-4 mA^2+s)^4);
  N[result, 600]
];

x41[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((J (7+J) (60 mA^2+(-23+J (7+J)) s))/(20 s^4 (-4 mA^2+s)^2));
  N[result, 600]
];

x50[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((412876800 mA^10+25804800 (-20+J (7+J)) mA^8 s+645120 (400+J (7+J) (-48+J (7+J))) mA^6 s^2+8960 (-7200+J (7+J) (1656+J (7+J) (-80+J (7+J)))) mA^4 s^3+80 (100800+J (7+J) (-44640+J (7+J) (3892+J (7+J) (-112+J (7+J))))) mA^2 s^4+(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) s^5)/(403200 s^6 (-4 mA^2+s)^5));
  N[result, 600]
];

x51[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (-2880 mA^4-36 (-64+3 J (7+J)) mA^2 s-(540+J (7+J) (-53+J (7+J))) s^2))/(360 s^5 (-4 mA^2+s)^3);
  N[result, 600]
];

x60[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((-89181388800 mA^12+5573836800 (24+J (7+J)) mA^10 s+139345920 (-600+J (7+J) (-58+J (7+J))) mA^8 s^2+1935360 (14400+J (7+J) (2520+(-7+J) J (7+J) (14+J))) mA^6 s^3+17280 (-302400+J (7+J) (-91008+J (7+J) (6132+J (7+J) (-140+J (7+J))))) mA^4 s^4+108 (4838400+J (7+J) (-58+J (7+J)) (-46080+J (7+J) (4152+J (7+J) (-122+J (7+J))))) mA^2 s^5+(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) s^6)/(21772800 s^7 (-4 mA^2+s)^6));
  N[result, 600]
];

x61[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (-806400 mA^6-8064 (-91+2 J (7+J)) mA^4 s-168 (1428+J (7+J) (-74+J (7+J))) mA^2 s^2-(-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) s^3))/(10080 s^6 (-4 mA^2+s)^4);
  N[result, 600]
];

x70[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((24970788864000 mA^14+1560674304000 (-28+J (7+J)) mA^12 s+39016857600 (840+J (7+J) (-68+J (7+J))) mA^10 s^2+541900800 (-25200+J (7+J) (3564+J (7+J) (-116+J (7+J)))) mA^8 s^3+4838400 (705600+J (7+J) (-88+J (7+J)) (1836+J (7+J) (-80+J (7+J)))) mA^6 s^4+30240 (-16934400+J (7+J) (6312960+J (7+J) (-532176+J (7+J) (16828+J (7+J) (-220+J (7+J)))))) mA^4 s^5+140 (304819200+J (7+J) (-203627520+J (7+J) (25466688+J (7+J) (-1218960+J (7+J) (26668+J (7+J) (-268+J (7+J))))))) mA^2 s^6+(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) s^7)/(1524096000 s^8 (-4 mA^2+s)^7));
  N[result, 600]
]

x71[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (-103219200 mA^8-3225600 (-40+J (7+J)) mA^6 s-17920 (3528+J (7+J) (-187+2 J (7+J))) mA^4 s^2-80 (-186336+J (7+J) (16156+J (7+J) (-392+3 J (7+J)))) mA^2 s^3-(1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) s^4))/(403200 s^7 (-4 mA^2+s)^5);
  N[result, 600]
];

x73[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (-7200 mA^4-180 (-28+J (7+J)) mA^2 s-(-2+J) (9+J) (-53+J (7+J)) s^2))/(360 s^7 (-4 mA^2+s)^3);
  N[result, 600]
];

x80[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((-8789717680128000 mA^16+549357355008000 (32+J (7+J)) mA^14 s+13733933875200 (-1120+(-6+J) J (7+J) (13+J)) mA^12 s^2+190749081600 (40320+J (7+J) (4788+J (7+J) (-134+J (7+J)))) mA^10 s^3+1703116800 (-1411200+J (7+J) (-261360+J (7+J) (12124+J (7+J) (-196+J (7+J))))) mA^8 s^4+10644480 (45158400+J (7+J) (12775680+J (7+J) (-887216+J (7+J) (23548+(-13+J) J (7+J) (20+J))))) mA^6 s^5+49280 (-1219276800+(-6+J) J (7+J) (13+J) (6981120+J (7+J) (-605424+J (7+J) (19516+J (7+J) (-244+J (7+J)))))) mA^4 s^6+176 (24385536000+J (7+J) (19294848000+J (7+J) (-2717588160+J (7+J) (150465168+J (7+J) (-4033640+J (7+J) (55608+J (7+J) (-378+J (7+J)))))))) mA^2 s^7+(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) s^8)/(134120448000 s^9 (-4 mA^2+s)^8));
  N[result, 600]
]

x81[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (-39016857600 mA^10-278691840 (-199+3 J (7+J)) mA^8 s-1935360 (16776+J (7+J) (-562+5 J (7+J))) mA^6 s^2-69120 (-143928+J (7+J) (8190+J (7+J) (-161+J (7+J)))) mA^4 s^3-108 (15298560+J (7+J) (-1351248+J (7+J) (44884+J (7+J) (-620+3 J (7+J))))) mA^2 s^4-(-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) s^5))/(21772800 s^8 (-4 mA^2+s)^6);
  N[result, 600]
];

x83[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (-1128960 mA^6-28224 (-38+J (7+J)) mA^4 s-56 (6516+J (7+J) (-382+5 J (7+J))) mA^2 s^2-(-2+J) (9+J) (2564+J (7+J) (-108+J (7+J))) s^3))/(5040 s^8 (-4 mA^2+s)^4);
  N[result, 600]
];

x90[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((3797158037815296000 mA^18+237322377363456000 (-36+J (7+J)) mA^16 s+5933059434086400 (1440+J (7+J) (-88+J (7+J))) mA^14 s^2+82403603251200 (-60480+J (7+J) (6192+J (7+J) (-152+J (7+J)))) mA^12 s^3+735746457600 (2540160+J (7+J) (-395424+J (7+J) (15876+J (7+J) (-224+J (7+J))))) mA^10 s^4+4598415360 (-101606400+J (7+J) (23230080+J (7+J) (-1372176+J (7+J) (31388+J (7+J) (-300+J (7+J)))))) mA^8 s^5+21288960 (3657830400+J (7+J) (-1234414080+J (7+J) (102113856+J (7+J) (-3399264+J (7+J) (52588+J (7+J) (-376+J (7+J))))))) mA^6 s^6+76032 (-109734912000+J (7+J) (57411763200+J (7+J) (-6511881600+J (7+J) (299402208+J (7+J) (-6732000+J (7+J) (78148+J (7+J) (-448+J (7+J)))))))) mA^4 s^7+216 (2414168064000+J (7+J) (-2228726016000+J (7+J) (345508727040+J (7+J) (-21390550272+J (7+J) (663989328+J (7+J) (-11259712+J (7+J) (105560+J (7+J) (-512+J (7+J))))))))) mA^2 s^8+(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^9)/(14485008384000 s^10 (-4 mA^2+s)^9));
  N[result, 600]
];

(* from x90 *)
largeJ[x_?NumericQ] := Module[{s, mA},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(14485008384000 (4 mA^2-s)^9 s);
  N[result, 600]
];

(* norm vector (snd state resonance) *)

g0snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (2 s)/(-4 mA^4+s^2);
  N[result, 600]
];

x10snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((8 mA^2+(-2+J (7+J)) s)/(2 s^2 (-4 mA^2+s)));
  N[result, 600]
];

x20snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((-320 mA^4+20 (8+J (7+J)) mA^2 s+(-20+J (7+J) (-13+J (7+J))) s^2)/(20 s^3 (-4 mA^2+s)^2));
  N[result, 600]
];

x30snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((23040 mA^6+1440 (-12+J (7+J)) mA^4 s+36 (120+J (7+J) (-28+J (7+J))) mA^2 s^2+(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) s^3)/(360 s^4 (-4 mA^2+s)^3));
  N[result, 600]
];

x40snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (2580480 mA^8-161280 (16+J (7+J)) mA^6 s-4032 (-240+J (7+J) (-38+J (7+J))) mA^4 s^2-56 (2880+J (7+J) (972+J (7+J) (-62+J (7+J)))) mA^2 s^3+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) s^4)/(10080 s^5 (-4 mA^2+s)^4);
  N[result, 600]
];

x41snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((J (7+J) (60 mA^2+(-23+J (7+J)) s))/(20 s^4 (-4 mA^2+s)^2));
  N[result, 600]
];

x50snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((412876800 mA^10+25804800 (-20+J (7+J)) mA^8 s+645120 (400+J (7+J) (-48+J (7+J))) mA^6 s^2+8960 (-7200+J (7+J) (1656+J (7+J) (-80+J (7+J)))) mA^4 s^3+80 (100800+J (7+J) (-44640+J (7+J) (3892+J (7+J) (-112+J (7+J))))) mA^2 s^4+(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) s^5)/(403200 s^6 (-4 mA^2+s)^5));
  N[result, 600]
];

x51snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7+J) (-2880 mA^4-36 (-64+3 J (7+J)) mA^2 s-(540+J (7+J) (-53+J (7+J))) s^2))/(360 s^5 (-4 mA^2+s)^3);
  N[result, 600]
];

x60snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((-89181388800 mA^12+5573836800 (24+J (7+J)) mA^10 s+139345920 (-600+J (7+J) (-58+J (7+J))) mA^8 s^2+1935360 (14400+J (7+J) (2520+(-7+J) J (7+J) (14+J))) mA^6 s^3+17280 (-302400+J (7+J) (-91008+J (7+J) (6132+J (7+J) (-140+J (7+J))))) mA^4 s^4+108 (4838400+J (7+J) (-58+J (7+J)) (-46080+J (7+J) (4152+J (7+J) (-122+J (7+J))))) mA^2 s^5+(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) s^6)/(21772800 s^7 (-4 mA^2+s)^6));
  N[result, 600]
];

x61snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7+J) (-806400 mA^6-8064 (-91+2 J (7+J)) mA^4 s-168 (1428+J (7+J) (-74+J (7+J))) mA^2 s^2-(-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) s^3))/(10080 s^6 (-4 mA^2+s)^4);
  N[result, 600]
];

x70snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((24970788864000 mA^14+1560674304000 (-28+J (7+J)) mA^12 s+39016857600 (840+J (7+J) (-68+J (7+J))) mA^10 s^2+541900800 (-25200+J (7+J) (3564+J (7+J) (-116+J (7+J)))) mA^8 s^3+4838400 (705600+J (7+J) (-88+J (7+J)) (1836+J (7+J) (-80+J (7+J)))) mA^6 s^4+30240 (-16934400+J (7+J) (6312960+J (7+J) (-532176+J (7+J) (16828+J (7+J) (-220+J (7+J)))))) mA^4 s^5+140 (304819200+J (7+J) (-203627520+J (7+J) (25466688+J (7+J) (-1218960+J (7+J) (26668+J (7+J) (-268+J (7+J))))))) mA^2 s^6+(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) s^7)/(1524096000 s^8 (-4 mA^2+s)^7));
  N[result, 600]
];

x71snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7+J) (-103219200 mA^8-3225600 (-40+J (7+J)) mA^6 s-17920 (3528+J (7+J) (-187+2 J (7+J))) mA^4 s^2-80 (-186336+J (7+J) (16156+J (7+J) (-392+3 J (7+J)))) mA^2 s^3-(1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) s^4))/(403200 s^7 (-4 mA^2+s)^5);
  N[result, 600]
];

x73snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7+J) (-7200 mA^4-180 (-28+J (7+J)) mA^2 s-(-2+J) (9+J) (-53+J (7+J)) s^2))/(360 s^7 (-4 mA^2+s)^3);
  N[result, 600]
];

x80snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((-8789717680128000 mA^16+549357355008000 (32+J (7+J)) mA^14 s+13733933875200 (-1120+(-6+J) J (7+J) (13+J)) mA^12 s^2+190749081600 (40320+J (7+J) (4788+J (7+J) (-134+J (7+J)))) mA^10 s^3+1703116800 (-1411200+J (7+J) (-261360+J (7+J) (12124+J (7+J) (-196+J (7+J))))) mA^8 s^4+10644480 (45158400+J (7+J) (12775680+J (7+J) (-887216+J (7+J) (23548+(-13+J) J (7+J) (20+J))))) mA^6 s^5+49280 (-1219276800+(-6+J) J (7+J) (13+J) (6981120+J (7+J) (-605424+J (7+J) (19516+J (7+J) (-244+J (7+J)))))) mA^4 s^6+176 (24385536000+J (7+J) (19294848000+J (7+J) (-2717588160+J (7+J) (150465168+J (7+J) (-4033640+J (7+J) (55608+J (7+J) (-378+J (7+J)))))))) mA^2 s^7+(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) s^8)/(134120448000 s^9 (-4 mA^2+s)^8));
  N[result, 600]
];

x81snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7+J) (-39016857600 mA^10-278691840 (-199+3 J (7+J)) mA^8 s-1935360 (16776+J (7+J) (-562+5 J (7+J))) mA^6 s^2-69120 (-143928+J (7+J) (8190+J (7+J) (-161+J (7+J)))) mA^4 s^3-108 (15298560+J (7+J) (-1351248+J (7+J) (44884+J (7+J) (-620+3 J (7+J))))) mA^2 s^4-(-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) s^5))/(21772800 s^8 (-4 mA^2+s)^6);
  N[result, 600]
];

x83snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7+J) (-1128960 mA^6-28224 (-38+J (7+J)) mA^4 s-56 (6516+J (7+J) (-382+5 J (7+J))) mA^2 s^2-(-2+J) (9+J) (2564+J (7+J) (-108+J (7+J))) s^3))/(5040 s^8 (-4 mA^2+s)^4);
  N[result, 600]
];

x90snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((3797158037815296000 mA^18+237322377363456000 (-36+J (7+J)) mA^16 s+5933059434086400 (1440+J (7+J) (-88+J (7+J))) mA^14 s^2+82403603251200 (-60480+J (7+J) (6192+J (7+J) (-152+J (7+J)))) mA^12 s^3+735746457600 (2540160+J (7+J) (-395424+J (7+J) (15876+J (7+J) (-224+J (7+J))))) mA^10 s^4+4598415360 (-101606400+J (7+J) (23230080+J (7+J) (-1372176+J (7+J) (31388+J (7+J) (-300+J (7+J)))))) mA^8 s^5+21288960 (3657830400+J (7+J) (-1234414080+J (7+J) (102113856+J (7+J) (-3399264+J (7+J) (52588+J (7+J) (-376+J (7+J))))))) mA^6 s^6+76032 (-109734912000+J (7+J) (57411763200+J (7+J) (-6511881600+J (7+J) (299402208+J (7+J) (-6732000+J (7+J) (78148+J (7+J) (-448+J (7+J)))))))) mA^4 s^7+216 (2414168064000+J (7+J) (-2228726016000+J (7+J) (345508727040+J (7+J) (-21390550272+J (7+J) (663989328+J (7+J) (-11259712+J (7+J) (105560+J (7+J) (-512+J (7+J))))))))) mA^2 s^8+(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^9)/(14485008384000 s^10 (-4 mA^2+s)^9));
  N[result, 600]
];


fList = {g0, x10, x20, x30, x40, x41, x50, x51, x60, x61, x70, x71, x73, x80, x81, x83, x90};
largeJList = {0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, largeJ};

(* norm = {G0, N2, N4, X52, X62, X72}; *)
norm = {g0snd, x10snd, x20snd, x30snd, x40snd, x41snd, x50snd, x51snd, x60snd, x61snd, x70snd, x71snd, x73snd, x80snd, x81snd, x83snd, x90snd};

Jmax  = 70;
Jlist = Range[0, Jmax, 2];   (* J = 0, 2, 4, …, 40 — exact discrete constraints *)

(* obj  = {-1, 0, 0, 0, 0, 0};   objective: maximise -y1 = minimise y1 *)

obj = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

xLeft  = SetPrecision[0, 600];   (* physical domain left endpoint  — check includes [xLeft,  x_min] *)
xRight = SetPrecision[1, 600];   (* physical domain right endpoint — check includes [x_max, xRight] *)

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


If[Length[largeJList] != Length[fList],
  Print["ERROR: largeJList has ", Length[largeJList],
        " entries but fList has ", Length[fList], ". They must match."];
  Quit[2]
];



(* F(x, J) = sum_{k=1}^{n} fList[[k]][x, J] * yVec[[k]]
   Evaluated for each (midpoint, J) pair in the interval check below. *)
F[x_?NumericQ, J_?IntegerQ] := Sum[yVec[[k]] * fList[[k]][x, J], {k, Length[yVec]}];

(* large J limit check *)
X[x_?NumericQ] := Sum[yVec[[k]] * largeJList[[k]][x], {k, Length[yVec]}];

(* Safety: clamp midpoints near known singularity at x -> 1 (sp = 1/(1-x))
   and provide a robust numeric evaluator that returns $Failed on
   non-finite or exceptional results. *)
singularTol = 10^-12;   (* distance from x=1 to avoid sp singularity *)

safeF[x_?NumericQ, J_?IntegerQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  val = Quiet[Check[N[F[x0, J], 600], $Failed]];
  If[val === $Failed || !NumberQ[val] ||
     MemberQ[{Infinity, -Infinity, ComplexInfinity, Indeterminate, DirectedInfinity}, val],
    $Failed,
    val
  ]
];

(* Safe evaluator for large-J functional X[x] *)
safeX[x_?NumericQ] := Module[{x0 = x, val},
  If[Abs[1 - x0] < singularTol, x0 = 1 - singularTol];
  val = Quiet[Check[N[X[x0], 600], $Failed]];
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
      Table[SetPrecision[xStar + (k - nPts/2) * s, 600], {k, 0, nPts}]
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

(* Print["Original points          : ", Length[samplePoints]];
Print["New refined points       : ", Length[newPoints]];
Print["Combined before dedup    : ", Length[allPoints]];
Print["Combined after  dedup    : ", Length[dedupPoints]]; *)
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
Print[""];
Print["STATUS: Refinement needed. Re-run pmp2sdp + sdpb with updated ", outFile, "."];

Quit[1];   (* exit code 1 = adaptive loop should continue *)
