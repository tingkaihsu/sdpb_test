(* New head for purely numerical data *)
ClearAll[NumericalPositiveMatrixWithPrefactor];

(* Helper: recursively stringify all numeric leaves *)
toJsonNestedNumberArray[expr_, prec_] := expr /. n_?NumericQ :> toJsonNumber[n, prec];

(*
  Design: the "constant-function" approach.

  For each sample point xᵢ we create ONE PositiveMatrixWithPrefactor block whose
  polynomial matrix entries are CONSTANT (degree-0) polynomials:
      W⁰ⱼ(x) = f1(xᵢ),   W¹ⱼ(x) = f2(xᵢ)

  This gives n independent 1×1 PSD constraints solved simultaneously:
      a·f1(xᵢ) + b·f2(xᵢ) ≥ 0   for each i = 1..n

  SDPB manual (Section 3.1, p.4) JSON nesting for "polynomials":
    Level 1: [ ]   — column list          (m_j columns; here m_j=1)
    Level 2:  [ ]  — row within column    (m_j rows;    here m_j=1)
    Level 3:   [ ] — polynomial vector    [Q^0, Q^1]  (N+1=2 entries for N=1)
    Level 4:    [ "c0", ..., "cd" ]       coefficient list of each polynomial

  For degree-0 (constant) polynomials dⱼ=0, so each coefficient list has 1 element:
    Level 4: [ "f1(xᵢ)" ]   and   [ "f2(xᵢ)" ]

  CORRECT JSON output for block i:
    "polynomials": [              ← level 1 (len=1)
        [                         ← level 2 (len=1)
            [                     ← level 3 (len=2): polynomial vector
                ["f1(xᵢ)"],       ← level 4 (len=1): Q^0 coefficient list
                ["f2(xᵢ)"]        ← level 4 (len=1): Q^1 coefficient list
            ]
        ]
    ]

  In Mathematica this requires exactly THREE wrapping brace pairs before the
  coefficient lists:
      {{{  {f1(xᵢ)}, {f2(xᵢ)}  }}}
       ↑1   ↑2        ↑↑ depth-3 elements = degree-0 coefficient lists

  FOUR wrapping pairs (the previous bug):
      {{{{  {f1(xᵢ)}, {f2(xᵢ)}  }}}}
  produced depth 5, pushing coefficient lists one level too deep. The parser,
  while in the Json_Float_Parser state (expecting float strings at depth 4),
  received another '[' (array start), triggering:
      "Not implemented: function 'bool json_start_array()' in class:
       '17Json_Float_ParserIN2El8BigFloatEE'"
*)

toJsonObject[NumericalPositiveMatrixWithPrefactor[pmp_?AssociationQ], prec_, getSampleDataFn_:Function[<||>]] :=
Module[
  {sampleData, samplePoints, sampleScalings, functionValues, basisValues,
   prefactor, reducedPrefactor},

  sampleData = getSampleDataFn[NumericalPositiveMatrixWithPrefactor[pmp], prec];

  prefactor        = Lookup[pmp, "prefactor",        1];
  reducedPrefactor = Lookup[pmp, "reducedPrefactor", prefactor];

  samplePoints   = Lookup[sampleData, "samplePoints",  Lookup[pmp, "samplePoints",  Missing[]]];
  sampleScalings = Lookup[sampleData, "sampleScalings", Lookup[pmp, "sampleScalings", Missing[]]];
  functionValues = Lookup[pmp, "polynomials", Missing[]];
  basisValues    = Lookup[sampleData, "basisValues", Lookup[pmp, "basisValues", Missing[]]];

  DeleteMissing @ <|
    "prefactor"        -> toJsonDampedRational[prefactor,        prec],
    "reducedPrefactor" -> toJsonDampedRational[reducedPrefactor, prec],
    "samplePoints"     -> toJsonNumberArray[samplePoints,   prec],
    "sampleScalings"   -> toJsonNumberArray[sampleScalings, prec],
    (* functionValues has structure {{{ {f1(xᵢ)}, {f2(xᵢ)} }}} — 3 structural
       wrapper levels, with depth-3 entries being the degree-0 coefficient lists.
       toJsonNestedNumberArray stringifies all numeric leaves in-place. *)
    "polynomials"      -> If[MissingQ[functionValues], Missing[],
      toJsonNestedNumberArray[functionValues, prec]
    ],
    "basisValues"      -> If[MissingQ[basisValues], Missing[], toJsonNestedNumberArray[basisValues, prec]]
  |>
];

WritePmpJsonNumerical[
  file_,
  SDP[objective_, normalization_, positiveMatricesWithPrefactors_],
  prec_,
  getSampleDataFn_:Function[<||>]
] := exportJson[
  file,
  <|
    "objective"     -> toJsonNumberArray[objective,     prec],
    "normalization" -> toJsonNumberArray[normalization, prec],
    "PositiveMatrixWithPrefactorArray" ->
      Table[toJsonObject[pmp, prec, getSampleDataFn], {pmp, positiveMatricesWithPrefactors}]
  |>
];


<< "../SDPB.m";

(* ================================================================
   PROBLEM-SPECIFIC SECTION  ← edit here for your problem
   ----------------------------------------------------------------
   Functions fk[x, J] depend on BOTH x (continuous, discretised)
   AND J (discrete spin, constrained EXACTLY for J in Jlist).
   Sampling is done over x; J constraints are imposed exactly.

   fList collects the functions for the generalised Table loop.
   Each fList[[k]][x, J] must return a numeric value.

   Jmax and Jlist define the discrete spin sum (even spins only).

   extraTriplet is the J → ∞ limiting constraint vector:
     As J → ∞, divide fk[x,J] by the leading power of J (here J^4):
       f1[x,J] / J^4 → 0        (f1 is J-independent)
       f2[x,J] / J^4 → 0        (f2 ~ J^2, subdominant)
       f3[x,J] / J^4 → 2        (f3 ~ 2·J^4, leading term)
     So extraTriplet = {0, 0, 2}.
     This enforces  0·y1 + 0·y2 + 2·y3 ≥ 0  in the J→∞ limit.
   ================================================================ *)
m1 = N[1/4, 600];
J1 = 0;
J2 = 2;
mgap = N[166/100, 600];
mAval = N[1/1000, 600];

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
  mA = N[maVal, 600];
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

(* first state resonance *)

g0fst[x_?NumericQ ] := Module[{s, mA, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  result = (2 s)/(-4 mA^4+s^2);
  N[result, 600]
];

x10fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((8 mA^2 + (-2 + J (7 + J)) s)/(2 s^2 (-4 mA^2 + s)));
  N[result, 600]
];

x20fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((-320 mA^4+20 (8+J (7+J)) mA^2 s+(-20+J (7+J) (-13+J (7+J))) s^2)/(20 s^3 (-4 mA^2+s)^2));
  N[result, 600]
];

x30fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((23040 mA^6+1440 (-12+J (7+J)) mA^4 s+36 (120+J (7+J) (-28+J (7+J))) mA^2 s^2+(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) s^3)/(360 s^4 (-4 mA^2+s)^3));
  N[result, 600]
];

x40fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (2580480 mA^8-161280 (16+J (7+J)) mA^6 s-4032 (-240+J (7+J) (-38+J (7+J))) mA^4 s^2-56 (2880+J (7+J) (972+J (7+J) (-62+J (7+J)))) mA^2 s^3+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) s^4)/(10080 s^5 (-4 mA^2+s)^4);
  N[result, 600]
];

x41fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((J (7+J) (60 mA^2+(-23+J (7+J)) s))/(20 s^4 (-4 mA^2+s)^2));
  N[result, 600]
];

x50fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((412876800 mA^10+25804800 (-20+J (7+J)) mA^8 s+645120 (400+J (7+J) (-48+J (7+J))) mA^6 s^2+8960 (-7200+J (7+J) (1656+J (7+J) (-80+J (7+J)))) mA^4 s^3+80 (100800+J (7+J) (-44640+J (7+J) (3892+J (7+J) (-112+J (7+J))))) mA^2 s^4+(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) s^5)/(403200 s^6 (-4 mA^2+s)^5));
  N[result, 600]
];

x51fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (-2880 mA^4-36 (-64+3 J (7+J)) mA^2 s-(540+J (7+J) (-53+J (7+J))) s^2))/(360 s^5 (-4 mA^2+s)^3);
  N[result, 600]
];

x60fst[x_?NumericQ ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((-89181388800 mA^12+5573836800 (24+J (7+J)) mA^10 s+139345920 (-600+J (7+J) (-58+J (7+J))) mA^8 s^2+1935360 (14400+J (7+J) (2520+(-7+J) J (7+J) (14+J))) mA^6 s^3+17280 (-302400+J (7+J) (-91008+J (7+J) (6132+J (7+J) (-140+J (7+J))))) mA^4 s^4+108 (4838400+J (7+J) (-58+J (7+J)) (-46080+J (7+J) (4152+J (7+J) (-122+J (7+J))))) mA^2 s^5+(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) s^6)/(21772800 s^7 (-4 mA^2+s)^6));
  N[result, 600]
];

x61fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (-806400 mA^6-8064 (-91+2 J (7+J)) mA^4 s-168 (1428+J (7+J) (-74+J (7+J))) mA^2 s^2-(-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) s^3))/(10080 s^6 (-4 mA^2+s)^4);
  N[result, 600]
];

x70fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((24970788864000 mA^14+1560674304000 (-28+J (7+J)) mA^12 s+39016857600 (840+J (7+J) (-68+J (7+J))) mA^10 s^2+541900800 (-25200+J (7+J) (3564+J (7+J) (-116+J (7+J)))) mA^8 s^3+4838400 (705600+J (7+J) (-88+J (7+J)) (1836+J (7+J) (-80+J (7+J)))) mA^6 s^4+30240 (-16934400+J (7+J) (6312960+J (7+J) (-532176+J (7+J) (16828+J (7+J) (-220+J (7+J)))))) mA^4 s^5+140 (304819200+J (7+J) (-203627520+J (7+J) (25466688+J (7+J) (-1218960+J (7+J) (26668+J (7+J) (-268+J (7+J))))))) mA^2 s^6+(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) s^7)/(1524096000 s^8 (-4 mA^2+s)^7));
  N[result, 600]
];

x71fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (-103219200 mA^8-3225600 (-40+J (7+J)) mA^6 s-17920 (3528+J (7+J) (-187+2 J (7+J))) mA^4 s^2-80 (-186336+J (7+J) (16156+J (7+J) (-392+3 J (7+J)))) mA^2 s^3-(1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) s^4))/(403200 s^7 (-4 mA^2+s)^5);
  N[result, 600]
];

x73fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (-7200 mA^4-180 (-28+J (7+J)) mA^2 s-(-2+J) (9+J) (-53+J (7+J)) s^2))/(360 s^7 (-4 mA^2+s)^3);
  N[result, 600]
];

x80fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((-8789717680128000 mA^16+549357355008000 (32+J (7+J)) mA^14 s+13733933875200 (-1120+(-6+J) J (7+J) (13+J)) mA^12 s^2+190749081600 (40320+J (7+J) (4788+J (7+J) (-134+J (7+J)))) mA^10 s^3+1703116800 (-1411200+J (7+J) (-261360+J (7+J) (12124+J (7+J) (-196+J (7+J))))) mA^8 s^4+10644480 (45158400+J (7+J) (12775680+J (7+J) (-887216+J (7+J) (23548+(-13+J) J (7+J) (20+J))))) mA^6 s^5+49280 (-1219276800+(-6+J) J (7+J) (13+J) (6981120+J (7+J) (-605424+J (7+J) (19516+J (7+J) (-244+J (7+J)))))) mA^4 s^6+176 (24385536000+J (7+J) (19294848000+J (7+J) (-2717588160+J (7+J) (150465168+J (7+J) (-4033640+J (7+J) (55608+J (7+J) (-378+J (7+J)))))))) mA^2 s^7+(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) s^8)/(134120448000 s^9 (-4 mA^2+s)^8));
  N[result, 600]
];

x81fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (-39016857600 mA^10-278691840 (-199+3 J (7+J)) mA^8 s-1935360 (16776+J (7+J) (-562+5 J (7+J))) mA^6 s^2-69120 (-143928+J (7+J) (8190+J (7+J) (-161+J (7+J)))) mA^4 s^3-108 (15298560+J (7+J) (-1351248+J (7+J) (44884+J (7+J) (-620+3 J (7+J))))) mA^2 s^4-(-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) s^5))/(21772800 s^8 (-4 mA^2+s)^6);
  N[result, 600]
];

x83fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (-1128960 mA^6-28224 (-38+J (7+J)) mA^4 s-56 (6516+J (7+J) (-382+5 J (7+J))) mA^2 s^2-(-2+J) (9+J) (2564+J (7+J) (-108+J (7+J))) s^3))/(5040 s^8 (-4 mA^2+s)^4);
  N[result, 600]
];

x90fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[maVal, 600];   (* Note: original x90 uses maVal, not mAval *)
  J = J1;
  result = -((3797158037815296000 mA^18+237322377363456000 (-36+J (7+J)) mA^16 s+5933059434086400 (1440+J (7+J) (-88+J (7+J))) mA^14 s^2+82403603251200 (-60480+J (7+J) (6192+J (7+J) (-152+J (7+J)))) mA^12 s^3+735746457600 (2540160+J (7+J) (-395424+J (7+J) (15876+J (7+J) (-224+J (7+J))))) mA^10 s^4+4598415360 (-101606400+J (7+J) (23230080+J (7+J) (-1372176+J (7+J) (31388+J (7+J) (-300+J (7+J)))))) mA^8 s^5+21288960 (3657830400+J (7+J) (-1234414080+J (7+J) (102113856+J (7+J) (-3399264+J (7+J) (52588+J (7+J) (-376+J (7+J))))))) mA^6 s^6+76032 (-109734912000+J (7+J) (57411763200+J (7+J) (-6511881600+J (7+J) (299402208+J (7+J) (-6732000+J (7+J) (78148+J (7+J) (-448+J (7+J)))))))) mA^4 s^7+216 (2414168064000+J (7+J) (-2228726016000+J (7+J) (345508727040+J (7+J) (-21390550272+J (7+J) (663989328+J (7+J) (-11259712+J (7+J) (105560+J (7+J) (-512+J (7+J))))))))) mA^2 s^8+(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^9)/(14485008384000 s^10 (-4 mA^2+s)^9));
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
  mA = N[maVal, 600];
  J = J2;
  result = -((3797158037815296000 mA^18+237322377363456000 (-36+J (7+J)) mA^16 s+5933059434086400 (1440+J (7+J) (-88+J (7+J))) mA^14 s^2+82403603251200 (-60480+J (7+J) (6192+J (7+J) (-152+J (7+J)))) mA^12 s^3+735746457600 (2540160+J (7+J) (-395424+J (7+J) (15876+J (7+J) (-224+J (7+J))))) mA^10 s^4+4598415360 (-101606400+J (7+J) (23230080+J (7+J) (-1372176+J (7+J) (31388+J (7+J) (-300+J (7+J)))))) mA^8 s^5+21288960 (3657830400+J (7+J) (-1234414080+J (7+J) (102113856+J (7+J) (-3399264+J (7+J) (52588+J (7+J) (-376+J (7+J))))))) mA^6 s^6+76032 (-109734912000+J (7+J) (57411763200+J (7+J) (-6511881600+J (7+J) (299402208+J (7+J) (-6732000+J (7+J) (78148+J (7+J) (-448+J (7+J)))))))) mA^4 s^7+216 (2414168064000+J (7+J) (-2228726016000+J (7+J) (345508727040+J (7+J) (-21390550272+J (7+J) (663989328+J (7+J) (-11259712+J (7+J) (105560+J (7+J) (-512+J (7+J))))))))) mA^2 s^8+(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^9)/(14485008384000 s^10 (-4 mA^2+s)^9));
  N[result, 600]
];


fList = {g0, x10, x20, x30, x40, x41, x50, x51, x60, x61, x70, x71, x73, x80, x81, x83, x90};
largeJList = {0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, largeJ};
ResList = {g0fst, x10fst, x20fst, x30fst, x40fst, x41fst, x50fst, x51fst, x60fst, x61fst, x70fst, x71fst, x73fst, x80fst, x81fst, x83fst, x90fst};

(* norm = {G0, N2, N4, X52, X62, X72}; *)
norm = {g0snd, x10snd, x20snd, x30snd, x40snd, x41snd, x50snd, x51snd, x60snd, x61snd, x70snd, x71snd, x73snd, x80snd, x81snd, x83snd, x90snd};

Jmax  = 70;
Jlist = Range[0, Jmax, 2];   (* J = 0, 2, 4, …, 40 — exact discrete constraints *)

(* obj  = {-1, 0, 0, 0, 0, 0};   objective: maximise -y1 = minimise y1 *)

obj = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

(* ================================================================
   END OF PROBLEM-SPECIFIC SECTION
   ================================================================ *)


(* ----------------------------------------------------------------
   testNumericalSDP
   ----------------------------------------------------------------
   Read x sample points from spFile (one number per line; blank
   lines and "#"-comment lines are ignored).

   For each (xi, Jj) pair a REGULAR block is created enforcing:
     f1(xi,Jj)·y1 + f2(xi,Jj)·y2 + f3(xi,Jj)·y3 ≥ 0

   For each xi an EXTRA block is created from extraTriplet,
   enforcing the J→∞ limit constraint:
     extraTriplet[[1]]·y1 + extraTriplet[[2]]·y2 + extraTriplet[[3]]·y3 ≥ 0

   Total blocks = Length[samplePoints] × (Length[Jlist] + 1).

   sampleScalings are Exp[-xi], the value of the prefactor
   DampedRational[1,{},1/E,x] at each sample point xi.
   ---------------------------------------------------------------- *)
testNumericalSDP[spFile_String, jsonFile_String, prec_:1000] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, polsRegular, polsExtra},

  (* --- Read and parse sampling_points.txt --- *)
  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " x sample points from ", spFile];
  (* Print["  x-points    : ", samplePoints];
  Print["  J-values    : ", Jlist, "  (", Length[Jlist], " spins, exact)"];
  Print["  extraTriplet: ", extraTriplet, "  (J\[Rule]\[Infinity] limit)"]; *)

  (* Scalings = prefactor DampedRational[1,{},1/E,x] evaluated at xi = e^{-xi} *)
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  (* --- Regular blocks: one per (xi, Jj) pair.
     Polynomials nesting:  {{ Table[{fk(xi,Jj)}, {k,3}] }}
       {{ ... }}  ← JSON levels 1 and 2  (column list / row, each size 1)
       Table[...] ← JSON level 3: polynomial vector, one entry per fList[[k]]
       {fk(xi,Jj)} ← JSON level 4: degree-0 coefficient list (1 element) --- *)

  pols = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[fList[[k]][samplePoints[[i]], Jlist[[j]]], prec]},
        {k, Length[fList]}
      ] }}
    |>],
    {i, Length[samplePoints]}, {j, Length[Jlist]}
  ];

  (* --- Extra blocks: one per xi, encoding the J→∞ limit constraint.
     extraTriplet = {0, 0, 2} is x-independent (J^4 leading coefficient),
     so every extra block carries the same polynomial values.
     Each block is still associated with a distinct xi for SDPB bookkeeping. --- *)
  polsLargeJ = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[largeJList[[k]][samplePoints[[i]]], prec]},
        {k, Length[largeJList]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  pols1Res = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[ResList[[k]][samplePoints[[i]]], prec]},
        {k, Length[ResList]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  (* Flatten polsRegular (2D Table → flat list) and append polsExtra (already flat) *)
  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Join[Flatten[pols], polsLargeJ, pols1Res]],
    prec
  ];
  Print["Wrote PMP JSON to ", jsonFile]
];


(* ----------------------------------------------------------------
   Command-line entry point
   ----------------------------------------------------------------
   USAGE:
     wolframscript -file g3_ExtremalEFT_2.m <sp_file> [output.json] [prec]

   ARGUMENTS:
     sp_file       required  path to sampling_points.txt (one x per line)
     output.json   optional  output path (default: numeric_pmp.json)
     prec          optional  decimal digit precision (default: 200)

   NOTE: $ScriptCommandLine = {scriptname, arg1, arg2, …} is populated
   identically by both wolframscript -file and math -script. Rest[] drops
   the script name leaving only user arguments. When loaded with << as a
   library, $ScriptCommandLine has length ≤ 1 and the block is skipped.
   ---------------------------------------------------------------- *)
Module[{myArgs, spFile, jsonFile, prec},

  myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

  If[Length[myArgs] >= 1,
    spFile   = myArgs[[1]];
    jsonFile = If[Length[myArgs] >= 2, myArgs[[2]], "n_pmp.json"];
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 600];

    Print["=== g3_ExtremalEFT_2.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    Null
  ]
];
