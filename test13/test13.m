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
m1 = N[3/4, 600];
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

x91[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(1524096000 s^9 (-4 mA^2+s)^7)J (7+J) (-9364045824000 mA^12-39016857600 (-416+7 J (7+J)) mA^10 s-3251404800 (3684+J (7+J) (-131+J (7+J))) mA^8 s^2-4838400 (-1005408+J (7+J) (57372+J (7+J) (-952+5 J (7+J)))) mA^6 s^3-120960 (9648000+J (7+J) (-798456+J (7+J) (21868+J (7+J) (-250+J (7+J))))) mA^4 s^4-420 (-387348480+J (7+J) (44625024+J (7+J) (-1824768+J (7+J) (34588+J (7+J) (-304+J (7+J)))))) mA^2 s^5-(11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))) s^6);
  N[result, 600]
];

x93[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(201600 s^9 (-4 mA^2+s)^5)J (7+J) (-361267200 mA^8-9031680 (-48+J (7+J)) mA^6 s-80640 (-6+J) (13+J) (-32+J (7+J)) mA^4 s^2-80 (-545760+J (7+J) (38892+J (7+J) (-784+5 J (7+J)))) mA^2 s^3-(-2+J) (9+J) (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))) s^4);
  N[result, 600]
];

x100[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -1/(1883051089920000 s^11 (-4 mA^2+s)^10)(-1974522179663953920000 mA^20+123407636228997120000 (40+J (7+J)) mA^18 s+3085190905724928000 (-1800+(-7+J) J (7+J) (14+J)) mA^16 s^2+42849873690624000 (86400+J (7+J) (7776+(-10+J) J (7+J) (17+J))) mA^14 s^3+382588157952000 (-4233600+J (7+J) (-568800+J (7+J) (20132+J (7+J) (-252+J (7+J))))) mA^12 s^4+2391175987200 (203212800+J (7+J) (39047040+J (7+J) (-2007216+J (7+J) (40348+J (7+J) (-340+J (7+J)))))) mA^10 s^5+11070259200 (-9144576000+J (7+J) (-2488838400+J (7+J) (176211360+J (7+J) (-5094216+J (7+J) (68788+J (7+J) (-430+J (7+J))))))) mA^8 s^6+39536640 (365783040000+J (7+J) (143820748800+J (7+J) (-13659851520+J (7+J) (537350688+J (7+J) (-10413160+J (7+J) (104468+J (7+J) (-518+J (7+J)))))))) mA^6 s^7+112320 (-12070840320000+J (7+J) (-7280961177600+J (7+J) (918554307840+J (7+J) (-47737944576+J (7+J) (1256405328+J (7+J) (-18136736+J (7+J) (144984+J (7+J) (-600+J (7+J))))))))) mA^4 s^8+260 (289700167680000+J (7+J) (308487979008000+J (7+J) (-51779666534400+J (7+J) (3504522827520+J (7+J) (-121704682752+J (7+J) (2396773008+J (7+J) (-27755072+J (7+J) (186600+J (7+J) (-672+J (7+J)))))))))) mA^2 s^9+(-1883051089920000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J))))))) s^10);
  N[result, 600]
];

x101[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(134120448000 s^10 (-4 mA^2+s)^8)J (7+J) (-4944216195072000 mA^14-27467867750400 (-347+4 J (7+J)) mA^12 s-190749081600 (41940+J (7+J) (-1046+7 J (7+J))) mA^10 s^2-3406233600 (-1119240+J (7+J) (45752+J (7+J) (-658+3 J (7+J)))) mA^8 s^3-53222400 (21139200+J (7+J) (-1275184+J (7+J) (29820+J (7+J) (-292+J (7+J))))) mA^6 s^4-98560 (-2123884800+J (7+J) (180272880+J (7+J) (-6162732+J (7+J) (98156+J (7+J) (-725+2 J (7+J)))))) mA^4 s^5-176 (134118374400+J (7+J) (-15741351360+J (7+J) (749269584+J (7+J) (-17497640+J (7+J) (211904+J (7+J) (-1274+3 J (7+J))))))) mA^2 s^6-(-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))) s^7);
  N[result, 600]
];

x103[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(10886400 s^10 (-4 mA^2+s)^6)J (7+J) (-133772083200 mA^10-3344302080 (-58+J (7+J)) mA^8 s-967680 (119088+J (7+J) (-4366+35 J (7+J))) mA^6 s^2-34560 (-1022544+J (7+J) (60018+5 J (7+J) (-203+J (7+J)))) mA^4 s^3-540 (10509696+J (7+J) (-885648+J (7+J) (24108+(-13+J) J (7+J) (20+J)))) mA^2 s^4-(-2+J) (9+J) (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))) s^5);
  N[result, 600]
];

x104[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(21772800 s^10 (-4 mA^2+s)^6)J (7+J) (-234101145600 mA^10-5852528640 (-58+J (7+J)) mA^8 s-54190080 (3708+J (7+J) (-134+J (7+J))) mA^6 s^2-17280 (-3527712+J (7+J) (200172+J (7+J) (-3080+13 J (7+J)))) mA^4 s^3-540 (17684352+J (7+J) (-1403760+J (7+J) (33852+J (7+J) (-308+J (7+J))))) mA^2 s^4-(-2+J) (9+J) (35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))) s^5);
  N[result, 600]
];

x110[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((1216305662672995614720000 mA^22+76019103917062225920000 (-4+J) (11+J) mA^20 s+1900477597926555648000 (2200+J (7+J) (-108+J (7+J))) mA^18 s^2+26395522193424384000 (-118800+J (7+J) (9540+J (7+J) (-188+J (7+J)))) mA^16 s^3+235674305298432000 (6652800+J (7+J) (-786528+J (7+J) (24892+J (7+J) (-280+J (7+J))))) mA^14 s^4+1472964408115200 (-372556800+J (7+J) (61799040+J (7+J) (-2812496+J (7+J) (50428+J (7+J) (-380+J (7+J)))))) mA^12 s^5+6819279667200 (20118067200+J (7+J) (-4597378560+J (7+J) (284601024+J (7+J) (-7273008+J (7+J) (87148+J (7+J) (-484+J (7+J))))))) mA^10 s^6+24354570240 (-1005903360000+J (7+J) (318039436800+J (7+J) (-25994646720+J (7+J) (893945808+J (7+J) (-15228320+J (7+J) (134568+(-21+J) J (7+J) (28+J))))))) mA^8 s^7+69189120 (44259747840000+J (7+J) (-19937187072000+J (7+J) (2120621241600+J (7+J) (-95024805120+J (7+J) (2172763408+J (7+J) (-27329920+J (7+J) (190568+J (7+J) (-688+J (7+J))))))))) mA^6 s^8+160160 (-1593350922240000+J (7+J) (1094831786188800+J (7+J) (-150983531781120+J (7+J) (8660220841728+J (7+J) (-257396458176+J (7+J) (4355540496+J (7+J) (-43413344+J (7+J) (251400+J (7+J) (-780+J (7+J)))))))))) mA^4 s^9+308 (-4+J) (11+J) (-941525544960000+(-10+J) (-8+J) (-6+J) (-2+J) J (7+J) (9+J) (13+J) (15+J) (17+J) (39263040+J (7+J) (-3046896+J (7+J) (60328+J (7+J) (-430+J (7+J)))))) mA^2 s^10+(-289989867847680000+(-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J))))))) s^11)/(289989867847680000 s^12 (-4 mA^2+s)^11));
  N[result, 600]
];

(* from x110 *)
largeJ[x_?NumericQ] := Module[{s, mA},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(289989867847680000 (4 mA^2-s)^11 s);
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
  mA = N[mAval, 600]; 
  J = J1;
  result = -((3797158037815296000 mA^18+237322377363456000 (-36+J (7+J)) mA^16 s+5933059434086400 (1440+J (7+J) (-88+J (7+J))) mA^14 s^2+82403603251200 (-60480+J (7+J) (6192+J (7+J) (-152+J (7+J)))) mA^12 s^3+735746457600 (2540160+J (7+J) (-395424+J (7+J) (15876+J (7+J) (-224+J (7+J))))) mA^10 s^4+4598415360 (-101606400+J (7+J) (23230080+J (7+J) (-1372176+J (7+J) (31388+J (7+J) (-300+J (7+J)))))) mA^8 s^5+21288960 (3657830400+J (7+J) (-1234414080+J (7+J) (102113856+J (7+J) (-3399264+J (7+J) (52588+J (7+J) (-376+J (7+J))))))) mA^6 s^6+76032 (-109734912000+J (7+J) (57411763200+J (7+J) (-6511881600+J (7+J) (299402208+J (7+J) (-6732000+J (7+J) (78148+J (7+J) (-448+J (7+J)))))))) mA^4 s^7+216 (2414168064000+J (7+J) (-2228726016000+J (7+J) (345508727040+J (7+J) (-21390550272+J (7+J) (663989328+J (7+J) (-11259712+J (7+J) (105560+J (7+J) (-512+J (7+J))))))))) mA^2 s^8+(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^9)/(14485008384000 s^10 (-4 mA^2+s)^9));
  N[result, 600]
];

x91fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(1524096000 s^9 (-4 mA^2+s)^7)J (7+J) (-9364045824000 mA^12-39016857600 (-416+7 J (7+J)) mA^10 s-3251404800 (3684+J (7+J) (-131+J (7+J))) mA^8 s^2-4838400 (-1005408+J (7+J) (57372+J (7+J) (-952+5 J (7+J)))) mA^6 s^3-120960 (9648000+J (7+J) (-798456+J (7+J) (21868+J (7+J) (-250+J (7+J))))) mA^4 s^4-420 (-387348480+J (7+J) (44625024+J (7+J) (-1824768+J (7+J) (34588+J (7+J) (-304+J (7+J)))))) mA^2 s^5-(11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))) s^6);
  N[result, 600]
];

x93fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(201600 s^9 (-4 mA^2+s)^5)J (7+J) (-361267200 mA^8-9031680 (-48+J (7+J)) mA^6 s-80640 (-6+J) (13+J) (-32+J (7+J)) mA^4 s^2-80 (-545760+J (7+J) (38892+J (7+J) (-784+5 J (7+J)))) mA^2 s^3-(-2+J) (9+J) (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))) s^4);
  N[result, 600]
];

x100fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -1/(1883051089920000 s^11 (-4 mA^2+s)^10)(-1974522179663953920000 mA^20+123407636228997120000 (40+J (7+J)) mA^18 s+3085190905724928000 (-1800+(-7+J) J (7+J) (14+J)) mA^16 s^2+42849873690624000 (86400+J (7+J) (7776+(-10+J) J (7+J) (17+J))) mA^14 s^3+382588157952000 (-4233600+J (7+J) (-568800+J (7+J) (20132+J (7+J) (-252+J (7+J))))) mA^12 s^4+2391175987200 (203212800+J (7+J) (39047040+J (7+J) (-2007216+J (7+J) (40348+J (7+J) (-340+J (7+J)))))) mA^10 s^5+11070259200 (-9144576000+J (7+J) (-2488838400+J (7+J) (176211360+J (7+J) (-5094216+J (7+J) (68788+J (7+J) (-430+J (7+J))))))) mA^8 s^6+39536640 (365783040000+J (7+J) (143820748800+J (7+J) (-13659851520+J (7+J) (537350688+J (7+J) (-10413160+J (7+J) (104468+J (7+J) (-518+J (7+J)))))))) mA^6 s^7+112320 (-12070840320000+J (7+J) (-7280961177600+J (7+J) (918554307840+J (7+J) (-47737944576+J (7+J) (1256405328+J (7+J) (-18136736+J (7+J) (144984+J (7+J) (-600+J (7+J))))))))) mA^4 s^8+260 (289700167680000+J (7+J) (308487979008000+J (7+J) (-51779666534400+J (7+J) (3504522827520+J (7+J) (-121704682752+J (7+J) (2396773008+J (7+J) (-27755072+J (7+J) (186600+J (7+J) (-672+J (7+J)))))))))) mA^2 s^9+(-1883051089920000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J))))))) s^10);
  N[result, 600]
];

x101fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(134120448000 s^10 (-4 mA^2+s)^8)J (7+J) (-4944216195072000 mA^14-27467867750400 (-347+4 J (7+J)) mA^12 s-190749081600 (41940+J (7+J) (-1046+7 J (7+J))) mA^10 s^2-3406233600 (-1119240+J (7+J) (45752+J (7+J) (-658+3 J (7+J)))) mA^8 s^3-53222400 (21139200+J (7+J) (-1275184+J (7+J) (29820+J (7+J) (-292+J (7+J))))) mA^6 s^4-98560 (-2123884800+J (7+J) (180272880+J (7+J) (-6162732+J (7+J) (98156+J (7+J) (-725+2 J (7+J)))))) mA^4 s^5-176 (134118374400+J (7+J) (-15741351360+J (7+J) (749269584+J (7+J) (-17497640+J (7+J) (211904+J (7+J) (-1274+3 J (7+J))))))) mA^2 s^6-(-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))) s^7);
  N[result, 600]
];

x103fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(10886400 s^10 (-4 mA^2+s)^6)J (7+J) (-133772083200 mA^10-3344302080 (-58+J (7+J)) mA^8 s-967680 (119088+J (7+J) (-4366+35 J (7+J))) mA^6 s^2-34560 (-1022544+J (7+J) (60018+5 J (7+J) (-203+J (7+J)))) mA^4 s^3-540 (10509696+J (7+J) (-885648+J (7+J) (24108+(-13+J) J (7+J) (20+J)))) mA^2 s^4-(-2+J) (9+J) (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))) s^5);
  N[result, 600]
];

x104fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(21772800 s^10 (-4 mA^2+s)^6)J (7+J) (-234101145600 mA^10-5852528640 (-58+J (7+J)) mA^8 s-54190080 (3708+J (7+J) (-134+J (7+J))) mA^6 s^2-17280 (-3527712+J (7+J) (200172+J (7+J) (-3080+13 J (7+J)))) mA^4 s^3-540 (17684352+J (7+J) (-1403760+J (7+J) (33852+J (7+J) (-308+J (7+J))))) mA^2 s^4-(-2+J) (9+J) (35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))) s^5);
  N[result, 600]
];

x110fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((1216305662672995614720000 mA^22+76019103917062225920000 (-4+J) (11+J) mA^20 s+1900477597926555648000 (2200+J (7+J) (-108+J (7+J))) mA^18 s^2+26395522193424384000 (-118800+J (7+J) (9540+J (7+J) (-188+J (7+J)))) mA^16 s^3+235674305298432000 (6652800+J (7+J) (-786528+J (7+J) (24892+J (7+J) (-280+J (7+J))))) mA^14 s^4+1472964408115200 (-372556800+J (7+J) (61799040+J (7+J) (-2812496+J (7+J) (50428+J (7+J) (-380+J (7+J)))))) mA^12 s^5+6819279667200 (20118067200+J (7+J) (-4597378560+J (7+J) (284601024+J (7+J) (-7273008+J (7+J) (87148+J (7+J) (-484+J (7+J))))))) mA^10 s^6+24354570240 (-1005903360000+J (7+J) (318039436800+J (7+J) (-25994646720+J (7+J) (893945808+J (7+J) (-15228320+J (7+J) (134568+(-21+J) J (7+J) (28+J))))))) mA^8 s^7+69189120 (44259747840000+J (7+J) (-19937187072000+J (7+J) (2120621241600+J (7+J) (-95024805120+J (7+J) (2172763408+J (7+J) (-27329920+J (7+J) (190568+J (7+J) (-688+J (7+J))))))))) mA^6 s^8+160160 (-1593350922240000+J (7+J) (1094831786188800+J (7+J) (-150983531781120+J (7+J) (8660220841728+J (7+J) (-257396458176+J (7+J) (4355540496+J (7+J) (-43413344+J (7+J) (251400+J (7+J) (-780+J (7+J)))))))))) mA^4 s^9+308 (-4+J) (11+J) (-941525544960000+(-10+J) (-8+J) (-6+J) (-2+J) J (7+J) (9+J) (13+J) (15+J) (17+J) (39263040+J (7+J) (-3046896+J (7+J) (60328+J (7+J) (-430+J (7+J)))))) mA^2 s^10+(-289989867847680000+(-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J))))))) s^11)/(289989867847680000 s^12 (-4 mA^2+s)^11));
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

x91snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(1524096000 s^9 (-4 mA^2 + s)^7) J (7 + J) (
     -9364045824000 mA^12
     - 39016857600 (-416 + 7 J (7 + J)) mA^10 s
     - 3251404800 (3684 + J (7 + J) (-131 + J (7 + J))) mA^8 s^2
     - 4838400 (-1005408 + J (7 + J) (57372 + J (7 + J) (-952 + 5 J (7 + J)))) mA^6 s^3
     - 120960 (9648000 + J (7 + J) (-798456 + J (7 + J) (21868 + J (7 + J) (-250 + J (7 + J))))) mA^4 s^4
     - 420 (-387348480 + J (7 + J) (44625024 + J (7 + J) (-1824768 + J (7 + J) (34588 + J (7 + J) (-304 + J (7 + J)))))) mA^2 s^5
     - (11405836800 + J (7 + J) (-1826254080 + J (7 + J) (107801568 + J (7 + J) (-3100260 + J (7 + J) (46228 + J (7 + J) (-343 + J (7 + J))))))) s^6
  );
  N[result, 600]
];

x93snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(201600 s^9 (-4 mA^2 + s)^5) J (7 + J) (
     -361267200 mA^8
     - 9031680 (-48 + J (7 + J)) mA^6 s
     - 80640 (-6 + J) (13 + J) (-32 + J (7 + J)) mA^4 s^2
     - 80 (-545760 + J (7 + J) (38892 + J (7 + J) (-784 + 5 J (7 + J)))) mA^2 s^3
     - (-2 + J) (9 + J) (-216000 + J (7 + J) (10752 + J (7 + J) (-182 + J (7 + J)))) s^4
  );
  N[result, 600]
];

x100snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -1/(1883051089920000 s^11 (-4 mA^2 + s)^10) (
     -1974522179663953920000 mA^20
     + 123407636228997120000 (40 + J (7 + J)) mA^18 s
     + 3085190905724928000 (-1800 + (-7 + J) J (7 + J) (14 + J)) mA^16 s^2
     + 42849873690624000 (86400 + J (7 + J) (7776 + (-10 + J) J (7 + J) (17 + J))) mA^14 s^3
     + 382588157952000 (-4233600 + J (7 + J) (-568800 + J (7 + J) (20132 + J (7 + J) (-252 + J (7 + J))))) mA^12 s^4
     + 2391175987200 (203212800 + J (7 + J) (39047040 + J (7 + J) (-2007216 + J (7 + J) (40348 + J (7 + J) (-340 + J (7 + J)))))) mA^10 s^5
     + 11070259200 (-9144576000 + J (7 + J) (-2488838400 + J (7 + J) (176211360 + J (7 + J) (-5094216 + J (7 + J) (68788 + J (7 + J) (-430 + J (7 + J))))))) mA^8 s^6
     + 39536640 (365783040000 + J (7 + J) (143820748800 + J (7 + J) (-13659851520 + J (7 + J) (537350688 + J (7 + J) (-10413160 + J (7 + J) (104468 + J (7 + J) (-518 + J (7 + J)))))))) mA^6 s^7
     + 112320 (-12070840320000 + J (7 + J) (-7280961177600 + J (7 + J) (918554307840 + J (7 + J) (-47737944576 + J (7 + J) (1256405328 + J (7 + J) (-18136736 + J (7 + J) (144984 + J (7 + J) (-600 + J (7 + J))))))))) mA^4 s^8
     + 260 (289700167680000 + J (7 + J) (308487979008000 + J (7 + J) (-51779666534400 + J (7 + J) (3504522827520 + J (7 + J) (-121704682752 + J (7 + J) (2396773008 + J (7 + J) (-27755072 + J (7 + J) (186600 + J (7 + J) (-672 + J (7 + J)))))))))) mA^2 s^9
     + (-1883051089920000 + (-8 + J) (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (15 + J) (-797572800 + J (7 + J) (106776000 + J (7 + J) (-3946556 + J (7 + J) (60108 + J (7 + J) (-405 + J (7 + J))))))) s^10
  );
  N[result, 600]
];

x101snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(134120448000 s^10 (-4 mA^2 + s)^8) J (7 + J) (
     -4944216195072000 mA^14
     - 27467867750400 (-347 + 4 J (7 + J)) mA^12 s
     - 190749081600 (41940 + J (7 + J) (-1046 + 7 J (7 + J))) mA^10 s^2
     - 3406233600 (-1119240 + J (7 + J) (45752 + J (7 + J) (-658 + 3 J (7 + J)))) mA^8 s^3
     - 53222400 (21139200 + J (7 + J) (-1275184 + J (7 + J) (29820 + J (7 + J) (-292 + J (7 + J))))) mA^6 s^4
     - 98560 (-2123884800 + J (7 + J) (180272880 + J (7 + J) (-6162732 + J (7 + J) (98156 + J (7 + J) (-725 + 2 J (7 + J)))))) mA^4 s^5
     - 176 (134118374400 + J (7 + J) (-15741351360 + J (7 + J) (749269584 + J (7 + J) (-17497640 + J (7 + J) (211904 + J (7 + J) (-1274 + 3 J (7 + J))))))) mA^2 s^6
     - (-1379752704000 + J (7 + J) (225934848000 + J (7 + J) (-14770082880 + J (7 + J) (486509168 + J (7 + J) (-8812960 + J (7 + J) (88928 + J (7 + J) (-468 + J (7 + J)))))))) s^7
  );
  N[result, 600]
];

x103snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(10886400 s^10 (-4 mA^2 + s)^6) J (7 + J) (
     -133772083200 mA^10
     - 3344302080 (-58 + J (7 + J)) mA^8 s
     - 967680 (119088 + J (7 + J) (-4366 + 35 J (7 + J))) mA^6 s^2
     - 34560 (-1022544 + J (7 + J) (60018 + 5 J (7 + J) (-203 + J (7 + J)))) mA^4 s^3
     - 540 (10509696 + J (7 + J) (-885648 + J (7 + J) (24108 + (-13 + J) J (7 + J) (20 + J)))) mA^2 s^4
     - (-2 + J) (9 + J) (21948480 + J (7 + J) (-1323000 + J (7 + J) (28702 + J (7 + J) (-277 + J (7 + J))))) s^5
  );
  N[result, 600]
];

x104snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(21772800 s^10 (-4 mA^2 + s)^6) J (7 + J) (
     -234101145600 mA^10
     - 5852528640 (-58 + J (7 + J)) mA^8 s
     - 54190080 (3708 + J (7 + J) (-134 + J (7 + J))) mA^6 s^2
     - 17280 (-3527712 + J (7 + J) (200172 + J (7 + J) (-3080 + 13 J (7 + J)))) mA^4 s^3
     - 540 (17684352 + J (7 + J) (-1403760 + J (7 + J) (33852 + J (7 + J) (-308 + J (7 + J))))) mA^2 s^4
     - (-2 + J) (9 + J) (35064000 + J (7 + J) (-1763640 + J (7 + J) (31942 + J (7 + J) (-277 + J (7 + J))))) s^5
  );
  N[result, 600]
];

x110snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -(
    (1216305662672995614720000 mA^22
      + 76019103917062225920000 (-4 + J) (11 + J) mA^20 s
      + 1900477597926555648000 (2200 + J (7 + J) (-108 + J (7 + J))) mA^18 s^2
      + 26395522193424384000 (-118800 + J (7 + J) (9540 + J (7 + J) (-188 + J (7 + J)))) mA^16 s^3
      + 235674305298432000 (6652800 + J (7 + J) (-786528 + J (7 + J) (24892 + J (7 + J) (-280 + J (7 + J))))) mA^14 s^4
      + 1472964408115200 (-372556800 + J (7 + J) (61799040 + J (7 + J) (-2812496 + J (7 + J) (50428 + J (7 + J) (-380 + J (7 + J)))))) mA^12 s^5
      + 6819279667200 (20118067200 + J (7 + J) (-4597378560 + J (7 + J) (284601024 + J (7 + J) (-7273008 + J (7 + J) (87148 + J (7 + J) (-484 + J (7 + J))))))) mA^10 s^6
      + 24354570240 (-1005903360000 + J (7 + J) (318039436800 + J (7 + J) (-25994646720 + J (7 + J) (893945808 + J (7 + J) (-15228320 + J (7 + J) (134568 + (-21 + J) J (7 + J) (28 + J))))))) mA^8 s^7
      + 69189120 (44259747840000 + J (7 + J) (-19937187072000 + J (7 + J) (2120621241600 + J (7 + J) (-95024805120 + J (7 + J) (2172763408 + J (7 + J) (-27329920 + J (7 + J) (190568 + J (7 + J) (-688 + J (7 + J))))))))) mA^6 s^8
      + 160160 (-1593350922240000 + J (7 + J) (1094831786188800 + J (7 + J) (-150983531781120 + J (7 + J) (8660220841728 + J (7 + J) (-257396458176 + J (7 + J) (4355540496 + J (7 + J) (-43413344 + J (7 + J) (251400 + J (7 + J) (-780 + J (7 + J)))))))))) mA^4 s^9
      + 308 (-4 + J) (11 + J) (-941525544960000 + (-10 + J) (-8 + J) (-6 + J) (-2 + J) J (7 + J) (9 + J) (13 + J) (15 + J) (17 + J) (39263040 + J (7 + J) (-3046896 + J (7 + J) (60328 + J (7 + J) (-430 + J (7 + J)))))) mA^2 s^10
      + (-289989867847680000 + (-10 + J) (-8 + J) (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (15 + J) (17 + J) (-833676480 + J (7 + J) (114378912 + J (7 + J) (-4230052 + J (7 + J) (63448 + J (7 + J) (-417 + J (7 + J))))))) s^11)
    /(289989867847680000 s^12 (-4 mA^2 + s)^11)
  );
  N[result, 600]
];

fList = {g0, x10, x20, x30, x40, x41, x50, x51, x60, x61, x70, x71, x73, x80, x81, x83, x90, x91, x93, x100, x101, x103, x104, x110};
largeJList = {0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, largeJ};
ResList = {g0fst, x10fst, x20fst, x30fst, x40fst, x41fst, x50fst, x51fst, x60fst, x61fst, x70fst, x71fst, x73fst, x80fst, x81fst, x83fst, x90fst, x91fst, x93fst, x100fst, x101fst, x103fst, x104fst, x110fst};

(* norm = {G0, N2, N4, X52, X62, X72}; *)
norm = {g0snd, x10snd, x20snd, x30snd, x40snd, x41snd, x50snd, x51snd, x60snd, x61snd, x70snd, x71snd, x73snd, x80snd, x81snd, x83snd, x90snd, x91snd, x93snd, x100snd, x101snd, x103snd, x104snd, x110snd};

Jmax  = 70;
Jlist = Range[0, Jmax, 2];   (* J = 0, 2, 4, …, 40 — exact discrete constraints *)

(* obj  = {-1, 0, 0, 0, 0, 0};   objective: maximise -y1 = minimise y1 *)

obj = {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

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
