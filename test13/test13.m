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
m1 = N[1/2, 600];
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
  result = 2/(-2 mA^2+s);
  N[result, 600]
];

x10[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (32 mA^4+2 (-1+J) (8+J) mA^2 s-(-2+J (7+J)) s^2)/(2 s^2 (-4 mA^2+s)^2);
  N[result, 600]
];

x20[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (-1280 mA^6+960 mA^4 s+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 s^2-(-20+J (7+J) (-13+J (7+J))) s^3)/(20 s^3 (-4 mA^2+s)^3);
  N[result, 600]
];

x30[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (92160 mA^8-92160 mA^6 s+34560 mA^4 s^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 s^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) s^4)/(360 s^4 (-4 mA^2+s)^4);
  N[result, 600]
];

x40[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(10080 s^5 (-4 mA^2+s)^5)(-10321920 mA^10+12902400 mA^8 s-6451200 mA^6 s^2+1612800 mA^4 s^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 s^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) s^5);
  N[result, 600]
];

x41[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (7+J) (11520 mA^8-11520 mA^6 s-4 (-936+J (7+J) (-26+J (7+J))) mA^4 s^2+2 (-216+J (7+J) (-26+J (7+J))) mA^2 s^3-9 (-23+J (7+J)) s^4))/(180 s^4 (-4 mA^2+s)^5);
  N[result, 600]
];

x50[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(403200 s^6 (-4 mA^2+s)^6)(1651507200 mA^12-2477260800 mA^10 s+1548288000 mA^8 s^2-516096000 mA^6 s^3+96768000 mA^4 s^4+2 (-4838400+(-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J)) mA^2 s^5-(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) s^6);
  N[result, 600]
];

x51[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(2520 s^5 (-4 mA^2+s)^6)J (7+J) (-645120 mA^10+806400 mA^8 s-403200 mA^6 s^2-2 (-54720+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^4 s^3+(-16920+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^2 s^4-7 (540+J (7+J) (-53+J (7+J))) s^5);
  N[result, 600]
];

x60[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(21772800 s^7 (-4 mA^2+s)^7)(-356725555200 mA^14+624269721600 mA^12 s-468202291200 mA^10 s^2+195084288000 mA^8 s^3-48771072000 mA^6 s^4+7315660800 mA^4 s^5+2 (-304819200+(-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J)) mA^2 s^6-(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) s^7);
  N[result, 600]
];

x61[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(100800 s^6 (-4 mA^2+s)^7)J (7+J) (103219200 mA^12-154828800 mA^10 s+96768000 mA^8 s^2-32256000 mA^6 s^3-2 (-2833920+J (7+J) (-44976+J (7+J) (3388+J (7+J) (-100+J (7+J))))) mA^4 s^4+(1+J) (6+J) (-69120+J (7+J) (4024+J (7+J) (-106+J (7+J)))) mA^2 s^5-10 (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) s^6);
  N[result, 600]
];

x70[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(1524096000 s^8 (-4 mA^2+s)^8)(99883155456000 mA^16-199766310912000 mA^14 s+174795522048000 mA^12 s^2-87397761024000 mA^10 s^3+27311800320000 mA^8 s^4-5462360064000 mA^6 s^5+682795008000 mA^4 s^6+2 (-24385536000+(-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J)) mA^2 s^7-(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) s^8);
  N[result, 600]
]

x71[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(10886400 s^7 (-4 mA^2+s)^8)J (7+J) (-44590694400 mA^14+78033715200 mA^12 s-58525286400 mA^10 s^2+24385536000 mA^8 s^3-6096384000 mA^6 s^4-4 (-240019200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^4 s^5+2 (-49507200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^2 s^6-27 (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) s^7);
  N[result, 600]
];

x73[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -1/(2520 s^5 (-4 mA^2+s)^8)(-2+J) J (7+J) (9+J) (3584 (-1+J) (8+J) mA^10+32 (-10+J) (-1+J) (8+J) (17+J) mA^8 s-32 (-1+J) (8+J) (-100+J (7+J)) mA^6 s^2+4 (-1+J) (8+J) (-230+3 J (7+J)) mA^4 s^3-2 (-1+J) (8+J) (-65+J (7+J)) mA^2 s^4+7 (-53+J (7+J)) s^5);
  N[result, 600]
];

x80[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(134120448000 s^9 (-4 mA^2+s)^9)(-35158870720512000 mA^18+79107459121152000 mA^16 s-79107459121152000 mA^14 s^2+46146017820672000 mA^12 s^3-17304756682752000 mA^10 s^4+4326189170688000 mA^8 s^5-721031528448000 mA^6 s^6+77253378048000 mA^4 s^7+2 (-2414168064000+(-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J)) mA^2 s^8-(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) s^9);
  N[result, 600]
]

x81[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(762048000 s^8 (-4 mA^2+s)^9)J (7+J) (12485394432000 mA^16-24970788864000 mA^14 s+21849440256000 mA^12 s^2-10924720128000 mA^10 s^3+3413975040000 mA^8 s^4-682795008000 mA^6 s^5-4 (-20447769600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^4 s^6+2 (-2158617600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^2 s^7-35 (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) s^8);
  N[result, 600]
];

x83[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(50400 s^6 (-4 mA^2+s)^9)(-2+J) J (7+J) (9+J) (286720 (-1+J) (8+J) mA^12-430080 (-1+J) (8+J) mA^10 s-16 (-1+J) (8+J) (-15480+J (7+J) (-74+J (7+J))) mA^8 s^2+16 (-1+J) (8+J) (-4280+J (7+J) (-74+J (7+J))) mA^6 s^3-6 (-1+J) (8+J) (-1480+J (7+J) (-74+J (7+J))) mA^4 s^4+(-1+J) (8+J) (-360+J (7+J) (-74+J (7+J))) mA^2 s^5-10 (2564+J (7+J) (-108+J (7+J))) s^6);
  N[result, 600]
];

x90[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(14485008384000 s^10 (-4 mA^2+s)^10)(15188632151261184000 mA^20-37971580378152960000 mA^18 s+42718027925422080000 mA^16 s^2-28478685283614720000 mA^14 s^3+12459424811581440000 mA^12 s^4-3737827443474432000 mA^10 s^5+778714050723840000 mA^8 s^6-111244864389120000 mA^6 s^7+10429206036480000 mA^4 s^8+2 (-289700167680000+(-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J)) mA^2 s^9-(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^10);
  N[result, 600]
];

(* x91[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(33530112000 s^9 (-4 mA^2+s)^10)J (7+J) (-2197429420032000 mA^18+4944216195072000 mA^16 s-4944216195072000 mA^14 s^2+2884126113792000 mA^12 s^3-1081547292672000 mA^10 s^4+270386823168000 mA^8 s^5-45064470528000 mA^6 s^6-2 (-2501346355200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^4 s^7+(-388949299200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^2 s^8-22 (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))) s^9);
  N[result, 600]
];

x93[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(5443200 s^7 (-4 mA^2+s)^10)(-2+J) J (7+J) (9+J) (-123863040 (-1+J) (8+J) mA^14+216760320 (-1+J) (8+J) mA^12 s-162570240 (-1+J) (8+J) mA^10 s^2-32 (-1+J) (8+J) (-2196000+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^8 s^3+32 (-1+J) (8+J) (-608400+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^6 s^4-12 (-1+J) (8+J) (-290880+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^4 s^5+2 (-1+J) (8+J) (-185040+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^2 s^6-27 (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))) s^7);
  N[result, 600]
];

x100[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(1883051089920000 s^11 (-4 mA^2+s)^11)(-7898088718655815680000 mA^22+21719743976303493120000 mA^20 s-27149679970379366400000 mA^18 s^2+20362259977784524800000 mA^16 s^3-10181129988892262400000 mA^14 s^4+3563395496112291840000 mA^12 s^5-890848874028072960000 mA^10 s^6+159080156076441600000 mA^8 s^7-19885019509555200000 mA^6 s^8+1657084959129600000 mA^4 s^9+2 (-41427123978240000+(-9+J) (-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J) (16+J)) mA^2 s^10+(1883051089920000-(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J))))))) s^11);
  N[result, 600]
];

x101[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(3621252096000 s^10 (-4 mA^2+s)^11)J (7+J) (949289509453824000 mA^20-2373223773634560000 mA^18 s+2669876745338880000 mA^16 s^2-1779917830225920000 mA^14 s^3+778714050723840000 mA^12 s^4-233614215217152000 mA^10 s^5+48669628170240000 mA^8 s^6-6952804024320000 mA^6 s^7-2 (-315451293696000+J (7+J) (-2977739366400+J (7+J) (314184925440+J (7+J) (-16618702464+J (7+J) (493173648+J (7+J) (-8546624+J (7+J) (85512+J (7+J) (-456+J (7+J))))))))) mA^4 s^8+(-25751126016000+J (7+J) (-2977739366400+J (7+J) (314184925440+J (7+J) (-16618702464+J (7+J) (493173648+J (7+J) (-8546624+J (7+J) (85512+J (7+J) (-456+J (7+J))))))))) mA^2 s^9-27 (-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))) s^10);
  N[result, 600]
];

x103[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(381024000 s^8 (-4 mA^2+s)^11)(-2+J) J (7+J) (9+J) (34681651200 (-1+J) (8+J) mA^16-69363302400 (-1+J) (8+J) mA^14 s+60692889600 (-1+J) (8+J) mA^12 s^2-30346444800 (-1+J) (8+J) mA^10 s^3-32 (-1+J) (8+J) (-290174400+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^8 s^4+32 (-1+J) (8+J) (-53092800+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^6 s^5-12 (-1+J) (8+J) (-13579200+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^4 s^6+2 (-1+J) (8+J) (-2289600+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^2 s^7-35 (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))) s^8);
  N[result, 600]
];

x104[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(21772800 s^7 (-4 mA^2+s)^11)(-2+J) J (7+J) (9+J) (-17694720 (-3+J) (-1+J) (8+J) (10+J) mA^14+30965760 (-3+J) (-1+J) (8+J) (10+J) mA^12 s+512 (-3+J) (-1+J) (8+J) (10+J) (-42720+J (7+J) (-104+J (7+J))) mA^10 s^2-640 (-3+J) (-1+J) (8+J) (10+J) (-12480+J (7+J) (-104+J (7+J))) mA^8 s^3+320 (-3+J) (-1+J) (8+J) (10+J) (-4920+J (7+J) (-104+J (7+J))) mA^6 s^4-80 (-3+J) (-1+J) (8+J) (10+J) (-1896+J (7+J) (-104+J (7+J))) mA^4 s^5+10 (-3+J) (-1+J) (8+J) (10+J) (-48+(-1+J) J) (8+J (15+J)) mA^2 s^6-(35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))) s^7);
  N[result, 600]
];

x110[x_?NumericQ, J_?IntegerQ] := Module[{s, mA, result},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = 1/(-4 mA^2+s)^12-1/(289989867847680000 (-4 mA^2+s)^12)(-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J))))))+(16777216 mA^24)/(s^12 (-4 mA^2+s)^12)-(50331648 mA^22)/(s^11 (-4 mA^2+s)^12)+(69206016 mA^20)/(s^10 (-4 mA^2+s)^12)-(57671680 mA^18)/(s^9 (-4 mA^2+s)^12)+(32440320 mA^16)/(s^8 (-4 mA^2+s)^12)-(12976128 mA^14)/(s^7 (-4 mA^2+s)^12)+(3784704 mA^12)/(s^6 (-4 mA^2+s)^12)-(811008 mA^10)/(s^5 (-4 mA^2+s)^12)+(126720 mA^8)/(s^4 (-4 mA^2+s)^12)-(14080 mA^6)/(s^3 (-4 mA^2+s)^12)+(1056 mA^4)/(s^2 (-4 mA^2+s)^12)+1/(144994933923840000 s (-4 mA^2+s)^12)(-6959756828344320000+(-10+J) (-9+J) (-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J) (16+J) (17+J)) mA^2;
  N[result, 600]
]; *)

(* from x80 *)
largeJ[x_?NumericQ] := Module[{s, mA},
  s = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((-2 mA^2+s)/(14485008384000 s (-4 mA^2+s)^10));
  N[result, 600]
];

(* first state resonance *)

(* heavy‑sum *)
g0fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 2/(-2 mA^2 + s);
  N[result, 600]
];

x10fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (32 mA^4+2 (-1+J) (8+J) mA^2 s-(-2+J (7+J)) s^2)/(2 s^2 (-4 mA^2+s)^2);
  N[result, 600]
];

x20fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (-1280 mA^6+960 mA^4 s+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 s^2-(-20+J (7+J) (-13+J (7+J))) s^3)/(20 s^3 (-4 mA^2+s)^3);
  N[result, 600]
];

x30fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (92160 mA^8-92160 mA^6 s+34560 mA^4 s^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 s^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) s^4)/(360 s^4 (-4 mA^2+s)^4);
  N[result, 600]
];

x40fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(10080 s^5 (-4 mA^2+s)^5)(-10321920 mA^10+12902400 mA^8 s-6451200 mA^6 s^2+1612800 mA^4 s^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 s^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) s^5);
  N[result, 600]
];

x41fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (7+J) (11520 mA^8-11520 mA^6 s-4 (-936+J (7+J) (-26+J (7+J))) mA^4 s^2+2 (-216+J (7+J) (-26+J (7+J))) mA^2 s^3-9 (-23+J (7+J)) s^4))/(180 s^4 (-4 mA^2+s)^5);
  N[result, 600]
];

x50fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(403200 s^6 (-4 mA^2+s)^6)(1651507200 mA^12-2477260800 mA^10 s+1548288000 mA^8 s^2-516096000 mA^6 s^3+96768000 mA^4 s^4+2 (-4838400+(-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J)) mA^2 s^5-(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) s^6);
  N[result, 600]
];

x51fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(2520 s^5 (-4 mA^2+s)^6)J (7+J) (-645120 mA^10+806400 mA^8 s-403200 mA^6 s^2-2 (-54720+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^4 s^3+(-16920+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^2 s^4-7 (540+J (7+J) (-53+J (7+J))) s^5);
  N[result, 600]
];

x60fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(21772800 s^7 (-4 mA^2+s)^7)(-356725555200 mA^14+624269721600 mA^12 s-468202291200 mA^10 s^2+195084288000 mA^8 s^3-48771072000 mA^6 s^4+7315660800 mA^4 s^5+2 (-304819200+(-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J)) mA^2 s^6-(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) s^7);
  N[result, 600]
];

x61fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(100800 s^6 (-4 mA^2+s)^7)J (7+J) (103219200 mA^12-154828800 mA^10 s+96768000 mA^8 s^2-32256000 mA^6 s^3-2 (-2833920+J (7+J) (-44976+J (7+J) (3388+J (7+J) (-100+J (7+J))))) mA^4 s^4+(1+J) (6+J) (-69120+J (7+J) (4024+J (7+J) (-106+J (7+J)))) mA^2 s^5-10 (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) s^6);
  N[result, 600]
];

x70fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(1524096000 s^8 (-4 mA^2+s)^8)(99883155456000 mA^16-199766310912000 mA^14 s+174795522048000 mA^12 s^2-87397761024000 mA^10 s^3+27311800320000 mA^8 s^4-5462360064000 mA^6 s^5+682795008000 mA^4 s^6+2 (-24385536000+(-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J)) mA^2 s^7-(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) s^8);
  N[result, 600]
];

x71fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(10886400 s^7 (-4 mA^2+s)^8)J (7+J) (-44590694400 mA^14+78033715200 mA^12 s-58525286400 mA^10 s^2+24385536000 mA^8 s^3-6096384000 mA^6 s^4-4 (-240019200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^4 s^5+2 (-49507200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^2 s^6-27 (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) s^7);
  N[result, 600]
];

x73fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -1/(2520 s^5 (-4 mA^2+s)^8)(-2+J) J (7+J) (9+J) (3584 (-1+J) (8+J) mA^10+32 (-10+J) (-1+J) (8+J) (17+J) mA^8 s-32 (-1+J) (8+J) (-100+J (7+J)) mA^6 s^2+4 (-1+J) (8+J) (-230+3 J (7+J)) mA^4 s^3-2 (-1+J) (8+J) (-65+J (7+J)) mA^2 s^4+7 (-53+J (7+J)) s^5);
  N[result, 600]
];

x80fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(134120448000 s^9 (-4 mA^2+s)^9)(-35158870720512000 mA^18+79107459121152000 mA^16 s-79107459121152000 mA^14 s^2+46146017820672000 mA^12 s^3-17304756682752000 mA^10 s^4+4326189170688000 mA^8 s^5-721031528448000 mA^6 s^6+77253378048000 mA^4 s^7+2 (-2414168064000+(-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J)) mA^2 s^8-(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) s^9);
  N[result, 600]
];

x81fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(762048000 s^8 (-4 mA^2+s)^9)J (7+J) (12485394432000 mA^16-24970788864000 mA^14 s+21849440256000 mA^12 s^2-10924720128000 mA^10 s^3+3413975040000 mA^8 s^4-682795008000 mA^6 s^5-4 (-20447769600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^4 s^6+2 (-2158617600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^2 s^7-35 (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) s^8);
  N[result, 600]
];

x83fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(50400 s^6 (-4 mA^2+s)^9)(-2+J) J (7+J) (9+J) (286720 (-1+J) (8+J) mA^12-430080 (-1+J) (8+J) mA^10 s-16 (-1+J) (8+J) (-15480+J (7+J) (-74+J (7+J))) mA^8 s^2+16 (-1+J) (8+J) (-4280+J (7+J) (-74+J (7+J))) mA^6 s^3-6 (-1+J) (8+J) (-1480+J (7+J) (-74+J (7+J))) mA^4 s^4+(-1+J) (8+J) (-360+J (7+J) (-74+J (7+J))) mA^2 s^5-10 (2564+J (7+J) (-108+J (7+J))) s^6);
  N[result, 600]
];

x90fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(14485008384000 s^10 (-4 mA^2+s)^10)(15188632151261184000 mA^20-37971580378152960000 mA^18 s+42718027925422080000 mA^16 s^2-28478685283614720000 mA^14 s^3+12459424811581440000 mA^12 s^4-3737827443474432000 mA^10 s^5+778714050723840000 mA^8 s^6-111244864389120000 mA^6 s^7+10429206036480000 mA^4 s^8+2 (-289700167680000+(-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J)) mA^2 s^9-(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) s^10);
  N[result, 600]
];

(* x91fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(33530112000 s^9 (-4 mA^2+s)^10)J (7+J) (-2197429420032000 mA^18+4944216195072000 mA^16 s-4944216195072000 mA^14 s^2+2884126113792000 mA^12 s^3-1081547292672000 mA^10 s^4+270386823168000 mA^8 s^5-45064470528000 mA^6 s^6-2 (-2501346355200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^4 s^7+(-388949299200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^2 s^8-22 (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))) s^9);
  N[result, 600]
];

x93fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(5443200 s^7 (-4 mA^2+s)^10)(-2+J) J (7+J) (9+J) (-123863040 (-1+J) (8+J) mA^14+216760320 (-1+J) (8+J) mA^12 s-162570240 (-1+J) (8+J) mA^10 s^2-32 (-1+J) (8+J) (-2196000+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^8 s^3+32 (-1+J) (8+J) (-608400+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^6 s^4-12 (-1+J) (8+J) (-290880+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^4 s^5+2 (-1+J) (8+J) (-185040+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^2 s^6-27 (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))) s^7);
  N[result, 600]
];

x100fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(1883051089920000 s^11 (-4 mA^2+s)^11)(-7898088718655815680000 mA^22+21719743976303493120000 mA^20 s-27149679970379366400000 mA^18 s^2+20362259977784524800000 mA^16 s^3-10181129988892262400000 mA^14 s^4+3563395496112291840000 mA^12 s^5-890848874028072960000 mA^10 s^6+159080156076441600000 mA^8 s^7-19885019509555200000 mA^6 s^8+1657084959129600000 mA^4 s^9+2 (-41427123978240000+(-9+J) (-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J) (16+J)) mA^2 s^10+(1883051089920000-(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J))))))) s^11);
  N[result, 600]
];

x101fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(3621252096000 s^10 (-4 mA^2+s)^11)J (7+J) (949289509453824000 mA^20-2373223773634560000 mA^18 s+2669876745338880000 mA^16 s^2-1779917830225920000 mA^14 s^3+778714050723840000 mA^12 s^4-233614215217152000 mA^10 s^5+48669628170240000 mA^8 s^6-6952804024320000 mA^6 s^7-2 (-315451293696000+J (7+J) (-2977739366400+J (7+J) (314184925440+J (7+J) (-16618702464+J (7+J) (493173648+J (7+J) (-8546624+J (7+J) (85512+J (7+J) (-456+J (7+J))))))))) mA^4 s^8+(-25751126016000+J (7+J) (-2977739366400+J (7+J) (314184925440+J (7+J) (-16618702464+J (7+J) (493173648+J (7+J) (-8546624+J (7+J) (85512+J (7+J) (-456+J (7+J))))))))) mA^2 s^9-27 (-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))) s^10);
  N[result, 600]
];

x103fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(381024000 s^8 (-4 mA^2+s)^11)(-2+J) J (7+J) (9+J) (34681651200 (-1+J) (8+J) mA^16-69363302400 (-1+J) (8+J) mA^14 s+60692889600 (-1+J) (8+J) mA^12 s^2-30346444800 (-1+J) (8+J) mA^10 s^3-32 (-1+J) (8+J) (-290174400+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^8 s^4+32 (-1+J) (8+J) (-53092800+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^6 s^5-12 (-1+J) (8+J) (-13579200+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^4 s^6+2 (-1+J) (8+J) (-2289600+J (7+J) (-528480+J (7+J) (16212+J (7+J) (-212+J (7+J))))) mA^2 s^7-35 (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))) s^8);
  N[result, 600]
];

x104fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(21772800 s^7 (-4 mA^2+s)^11)(-2+J) J (7+J) (9+J) (-17694720 (-3+J) (-1+J) (8+J) (10+J) mA^14+30965760 (-3+J) (-1+J) (8+J) (10+J) mA^12 s+512 (-3+J) (-1+J) (8+J) (10+J) (-42720+J (7+J) (-104+J (7+J))) mA^10 s^2-640 (-3+J) (-1+J) (8+J) (10+J) (-12480+J (7+J) (-104+J (7+J))) mA^8 s^3+320 (-3+J) (-1+J) (8+J) (10+J) (-4920+J (7+J) (-104+J (7+J))) mA^6 s^4-80 (-3+J) (-1+J) (8+J) (10+J) (-1896+J (7+J) (-104+J (7+J))) mA^4 s^5+10 (-3+J) (-1+J) (8+J) (10+J) (-48+(-1+J) J) (8+J (15+J)) mA^2 s^6-(35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))) s^7);
  N[result, 600]
];

x110fst[x_?NumericQ] := Module[{s, mA, J, result},
  s = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = 1/(-4 mA^2+s)^12-1/(289989867847680000 (-4 mA^2+s)^12)(-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J))))))+(16777216 mA^24)/(s^12 (-4 mA^2+s)^12)-(50331648 mA^22)/(s^11 (-4 mA^2+s)^12)+(69206016 mA^20)/(s^10 (-4 mA^2+s)^12)-(57671680 mA^18)/(s^9 (-4 mA^2+s)^12)+(32440320 mA^16)/(s^8 (-4 mA^2+s)^12)-(12976128 mA^14)/(s^7 (-4 mA^2+s)^12)+(3784704 mA^12)/(s^6 (-4 mA^2+s)^12)-(811008 mA^10)/(s^5 (-4 mA^2+s)^12)+(126720 mA^8)/(s^4 (-4 mA^2+s)^12)-(14080 mA^6)/(s^3 (-4 mA^2+s)^12)+(1056 mA^4)/(s^2 (-4 mA^2+s)^12)+1/(144994933923840000 s (-4 mA^2+s)^12)(-6959756828344320000+(-10+J) (-9+J) (-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J) (16+J) (17+J)) mA^2;
  N[result, 600]
]; *)

(* second state resonance *)
(* as norm vector *)
g0snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 2/(-2 mA^2 + s);
  N[result, 600]
];

x10snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (32 mA^4+2 (-1+J) (8+J) mA^2 s-(-2+J (7+J)) s^2)/(2 s^2 (-4 mA^2+s)^2);
  N[result, 600]
];

x20snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (-1280 mA^6 + 960 mA^4 s + 2 (-120 + (-1 + J) J (7 + J) (8 + J)) mA^2 s^2 - (-20 + J (7 + J) (-13 + J (7 + J))) s^3)/(20 s^3 (-4 mA^2 + s)^3);
  N[result, 600]
];

x30snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (92160 mA^8 - 92160 mA^6 s + 34560 mA^4 s^2 + 2 (-2880 + (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J)) mA^2 s^3 - (-12 + J (7 + J)) (30 + J (7 + J) (-23 + J (7 + J))) s^4)/(360 s^4 (-4 mA^2 + s)^4);
  N[result, 600]
];

x40snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(10080 s^5 (-4 mA^2 + s)^5) (-10321920 mA^10 + 12902400 mA^8 s - 6451200 mA^6 s^2 + 1612800 mA^4 s^3 + 2 (-100800 + (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J)) mA^2 s^4 + (10080 - (-2 + J) J (7 + J) (9 + J) (604 + J (7 + J) (-52 + J (7 + J)))) s^5);
  N[result, 600]
];

x41snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (7 + J) (11520 mA^8 - 11520 mA^6 s - 4 (-936 + J (7 + J) (-26 + J (7 + J))) mA^4 s^2 + 2 (-216 + J (7 + J) (-26 + J (7 + J))) mA^2 s^3 - 9 (-23 + J (7 + J)) s^4))/(180 s^4 (-4 mA^2 + s)^5);
  N[result, 600]
];

x50snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(403200 s^6 (-4 mA^2 + s)^6) (1651507200 mA^12 - 2477260800 mA^10 s + 1548288000 mA^8 s^2 - 516096000 mA^6 s^3 + 96768000 mA^4 s^4 + 2 (-4838400 + (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J)) mA^2 s^5 - (-403200 + (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (-34 + J (5 + J)) (-20 + J (9 + J))) s^6);
  N[result, 600]
];

x51snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(2520 s^5 (-4 mA^2 + s)^6) J (7 + J) (-645120 mA^10 + 806400 mA^8 s - 403200 mA^6 s^2 - 2 (-54720 + J (7 + J) (924 + J (7 + J) (-56 + J (7 + J)))) mA^4 s^3 + (-16920 + J (7 + J) (924 + J (7 + J) (-56 + J (7 + J)))) mA^2 s^4 - 7 (540 + J (7 + J) (-53 + J (7 + J))) s^5);
  N[result, 600]
];

x60snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(21772800 s^7 (-4 mA^2 + s)^7) (-356725555200 mA^14 + 624269721600 mA^12 s - 468202291200 mA^10 s^2 + 195084288000 mA^8 s^3 - 48771072000 mA^6 s^4 + 7315660800 mA^4 s^5 + 2 (-304819200 + (-5 + J) (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J) (12 + J)) mA^2 s^6 - (-21772800 + (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (-62 + J (7 + J)) (-48 + J (7 + J)) (-15 + J (7 + J))) s^7);
  N[result, 600]
];

x61snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(100800 s^6 (-4 mA^2 + s)^7) J (7 + J) (103219200 mA^12 - 154828800 mA^10 s + 96768000 mA^8 s^2 - 32256000 mA^6 s^3 - 2 (-2833920 + J (7 + J) (-44976 + J (7 + J) (3388 + J (7 + J) (-100 + J (7 + J))))) mA^4 s^4 + (1 + J) (6 + J) (-69120 + J (7 + J) (4024 + J (7 + J) (-106 + J (7 + J)))) mA^2 s^5 - 10 (-31032 + J (7 + J) (3024 + (-7 + J) J (7 + J) (14 + J))) s^6);
  N[result, 600]
];

x70snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(1524096000 s^8 (-4 mA^2 + s)^8) (99883155456000 mA^16 - 199766310912000 mA^14 s + 174795522048000 mA^12 s^2 - 87397761024000 mA^10 s^3 + 27311800320000 mA^8 s^4 - 5462360064000 mA^6 s^5 + 682795008000 mA^4 s^6 + 2 (-24385536000 + (-6 + J) (-5 + J) (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J) (12 + J) (13 + J)) mA^2 s^7 - (-1524096000 + (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (-50 + J (7 + J)) (960 + J (7 + J) (-83 + J (7 + J)))) s^8);
  N[result, 600]
];

x71snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(10886400 s^7 (-4 mA^2 + s)^8) J (7 + J) (-44590694400 mA^14 + 78033715200 mA^12 s - 58525286400 mA^10 s^2 + 24385536000 mA^8 s^3 - 6096384000 mA^6 s^4 - 4 (-240019200 + J (7 + J) (2888640 + J (7 + J) (-248256 + J (7 + J) (9388 + J (7 + J) (-160 + J (7 + J)))))) mA^4 s^5 + 2 (-49507200 + J (7 + J) (2888640 + J (7 + J) (-248256 + J (7 + J) (9388 + J (7 + J) (-160 + J (7 + J)))))) mA^2 s^6 - 27 (1578240 + J (7 + J) (-209056 + J (7 + J) (8988 + J (7 + J) (-160 + J (7 + J))))) s^7);
  N[result, 600]
];

x73snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -1/(2520 s^5 (-4 mA^2 + s)^8) (-2 + J) J (7 + J) (9 + J) (3584 (-1 + J) (8 + J) mA^10 + 32 (-10 + J) (-1 + J) (8 + J) (17 + J) mA^8 s - 32 (-1 + J) (8 + J) (-100 + J (7 + J)) mA^6 s^2 + 4 (-1 + J) (8 + J) (-230 + 3 J (7 + J)) mA^4 s^3 - 2 (-1 + J) (8 + J) (-65 + J (7 + J)) mA^2 s^4 + 7 (-53 + J (7 + J)) s^5);
  N[result, 600]
];

x80snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(134120448000 s^9 (-4 mA^2 + s)^9) (-35158870720512000 mA^18 + 79107459121152000 mA^16 s - 79107459121152000 mA^14 s^2 + 46146017820672000 mA^12 s^3 - 17304756682752000 mA^10 s^4 + 4326189170688000 mA^8 s^5 - 721031528448000 mA^6 s^6 + 77253378048000 mA^4 s^7 + 2 (-2414168064000 + (-7 + J) (-6 + J) (-5 + J) (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J) (12 + J) (13 + J) (14 + J)) mA^2 s^8 - (-134120448000 + (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (5001600 + J (7 + J) (-600160 + J (7 + J) (19516 + J (7 + J) (-240 + J (7 + J)))))) s^9);
  N[result, 600]
];

x81snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(762048000 s^8 (-4 mA^2 + s)^9) J (7 + J) (12485394432000 mA^16 - 24970788864000 mA^14 s + 21849440256000 mA^12 s^2 - 10924720128000 mA^10 s^3 + 3413975040000 mA^8 s^4 - 682795008000 mA^6 s^5 - 4 (-20447769600 + J (7 + J) (-236718720 + J (7 + J) (22252608 + J (7 + J) (-980520 + J (7 + J) (21868 + J (7 + J) (-238 + J (7 + J))))))) mA^4 s^6 + 2 (-2158617600 + J (7 + J) (-236718720 + J (7 + J) (22252608 + J (7 + J) (-980520 + J (7 + J) (21868 + J (7 + J) (-238 + J (7 + J))))))) mA^2 s^7 - 35 (-131466240 + J (7 + J) (17720496 + J (7 + J) (-915804 + J (7 + J) (21808 + J (7 + J) (-241 + J (7 + J)))))) s^8);
  N[result, 600]
];

x83snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(50400 s^6 (-4 mA^2 + s)^9) (-2 + J) J (7 + J) (9 + J) (286720 (-1 + J) (8 + J) mA^12 - 430080 (-1 + J) (8 + J) mA^10 s - 16 (-1 + J) (8 + J) (-15480 + J (7 + J) (-74 + J (7 + J))) mA^8 s^2 + 16 (-1 + J) (8 + J) (-4280 + J (7 + J) (-74 + J (7 + J))) mA^6 s^3 - 6 (-1 + J) (8 + J) (-1480 + J (7 + J) (-74 + J (7 + J))) mA^4 s^4 + (-1 + J) (8 + J) (-360 + J (7 + J) (-74 + J (7 + J))) mA^2 s^5 - 10 (2564 + J (7 + J) (-108 + J (7 + J))) s^6);
  N[result, 600]
];

x90snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(14485008384000 s^10 (-4 mA^2 + s)^10) (15188632151261184000 mA^20 - 37971580378152960000 mA^18 s + 42718027925422080000 mA^16 s^2 - 28478685283614720000 mA^14 s^3 + 12459424811581440000 mA^12 s^4 - 3737827443474432000 mA^10 s^5 + 778714050723840000 mA^8 s^6 - 111244864389120000 mA^6 s^7 + 10429206036480000 mA^4 s^8 + 2 (-289700167680000 + (-8 + J) (-7 + J) (-6 + J) (-5 + J) (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J) (12 + J) (13 + J) (14 + J) (15 + J)) mA^2 s^9 - (-14485008384000 + (-8 + J) (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (15 + J) (5277600 + J (7 + J) (-651672 + J (7 + J) (20980 + J (7 + J) (-250 + J (7 + J)))))) s^10);
  N[result, 600]
];

(* x91snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(33530112000 s^9 (-4 mA^2 + s)^10) J (7 + J) (-2197429420032000 mA^18 + 4944216195072000 mA^16 s - 4944216195072000 mA^14 s^2 + 2884126113792000 mA^12 s^3 - 1081547292672000 mA^10 s^4 + 270386823168000 mA^8 s^5 - 45064470528000 mA^6 s^6 - 2 (-2501346355200 + J (7 + J) (24088008960 + J (7 + J) (-2417474304 + J (7 + J) (118343568 + J (7 + J) (-3123584 + J (7 + J) (45192 + J (7 + J) (-336 + J (7 + J)))))))) mA^4 s^7 + (-388949299200 + J (7 + J) (24088008960 + J (7 + J) (-2417474304 + J (7 + J) (118343568 + J (7 + J) (-3123584 + J (7 + J) (45192 + J (7 + J) (-336 + J (7 + J)))))))) mA^2 s^8 - 22 (11405836800 + J (7 + J) (-1826254080 + J (7 + J) (107801568 + J (7 + J) (-3100260 + J (7 + J) (46228 + J (7 + J) (-343 + J (7 + J))))))) s^9);
  N[result, 600]
];

x93snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(5443200 s^7 (-4 mA^2 + s)^10) (-2 + J) J (7 + J) (9 + J) (-123863040 (-1 + J) (8 + J) mA^14 + 216760320 (-1 + J) (8 + J) mA^12 s - 162570240 (-1 + J) (8 + J) mA^10 s^2 - 32 (-1 + J) (8 + J) (-2196000 + J (7 + J) (5760 + J (7 + J) (-134 + J (7 + J)))) mA^8 s^3 + 32 (-1 + J) (8 + J) (-608400 + J (7 + J) (5760 + J (7 + J) (-134 + J (7 + J)))) mA^6 s^4 - 12 (-1 + J) (8 + J) (-290880 + J (7 + J) (5760 + J (7 + J) (-134 + J (7 + J)))) mA^4 s^5 + 2 (-1 + J) (8 + J) (-185040 + J (7 + J) (5760 + J (7 + J) (-134 + J (7 + J)))) mA^2 s^6 - 27 (-216000 + J (7 + J) (10752 + J (7 + J) (-182 + J (7 + J)))) s^7);
  N[result, 600]
];

x100snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(1883051089920000 s^11 (-4 mA^2 + s)^11) (-7898088718655815680000 mA^22 + 21719743976303493120000 mA^20 s - 27149679970379366400000 mA^18 s^2 + 20362259977784524800000 mA^16 s^3 - 10181129988892262400000 mA^14 s^4 + 3563395496112291840000 mA^12 s^5 - 890848874028072960000 mA^10 s^6 + 159080156076441600000 mA^8 s^7 - 19885019509555200000 mA^6 s^8 + 1657084959129600000 mA^4 s^9 + 2 (-41427123978240000 + (-9 + J) (-8 + J) (-7 + J) (-6 + J) (-5 + J) (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J) (12 + J) (13 + J) (14 + J) (15 + J) (16 + J)) mA^2 s^10 + (1883051089920000 - (-8 + J) (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (15 + J) (-797572800 + J (7 + J) (106776000 + J (7 + J) (-3946556 + J (7 + J) (60108 + J (7 + J) (-405 + J (7 + J))))))) s^11);
  N[result, 600]
];

x101snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(3621252096000 s^10 (-4 mA^2 + s)^11) J (7 + J) (949289509453824000 mA^20 - 2373223773634560000 mA^18 s + 2669876745338880000 mA^16 s^2 - 1779917830225920000 mA^14 s^3 + 778714050723840000 mA^12 s^4 - 233614215217152000 mA^10 s^5 + 48669628170240000 mA^8 s^6 - 6952804024320000 mA^6 s^7 - 2 (-315451293696000 + J (7 + J) (-2977739366400 + J (7 + J) (314184925440 + J (7 + J) (-16618702464 + J (7 + J) (493173648 + J (7 + J) (-8546624 + J (7 + J) (85512 + J (7 + J) (-456 + J (7 + J))))))))) mA^4 s^8 + (-25751126016000 + J (7 + J) (-2977739366400 + J (7 + J) (314184925440 + J (7 + J) (-16618702464 + J (7 + J) (493173648 + J (7 + J) (-8546624 + J (7 + J) (85512 + J (7 + J) (-456 + J (7 + J))))))))) mA^2 s^9 - 27 (-1379752704000 + J (7 + J) (225934848000 + J (7 + J) (-14770082880 + J (7 + J) (486509168 + J (7 + J) (-8812960 + J (7 + J) (88928 + J (7 + J) (-468 + J (7 + J)))))))) s^10);
  N[result, 600]
];

x103snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(381024000 s^8 (-4 mA^2 + s)^11) (-2 + J) J (7 + J) (9 + J) (34681651200 (-1 + J) (8 + J) mA^16 - 69363302400 (-1 + J) (8 + J) mA^14 s + 60692889600 (-1 + J) (8 + J) mA^12 s^2 - 30346444800 (-1 + J) (8 + J) mA^10 s^3 - 32 (-1 + J) (8 + J) (-290174400 + J (7 + J) (-528480 + J (7 + J) (16212 + J (7 + J) (-212 + J (7 + J))))) mA^8 s^4 + 32 (-1 + J) (8 + J) (-53092800 + J (7 + J) (-528480 + J (7 + J) (16212 + J (7 + J) (-212 + J (7 + J))))) mA^6 s^5 - 12 (-1 + J) (8 + J) (-13579200 + J (7 + J) (-528480 + J (7 + J) (16212 + J (7 + J) (-212 + J (7 + J))))) mA^4 s^6 + 2 (-1 + J) (8 + J) (-2289600 + J (7 + J) (-528480 + J (7 + J) (16212 + J (7 + J) (-212 + J (7 + J))))) mA^2 s^7 - 35 (21948480 + J (7 + J) (-1323000 + J (7 + J) (28702 + J (7 + J) (-277 + J (7 + J))))) s^8);
  N[result, 600]
];

x104snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(21772800 s^7 (-4 mA^2 + s)^11) (-2 + J) J (7 + J) (9 + J) (-17694720 (-3 + J) (-1 + J) (8 + J) (10 + J) mA^14 + 30965760 (-3 + J) (-1 + J) (8 + J) (10 + J) mA^12 s + 512 (-3 + J) (-1 + J) (8 + J) (10 + J) (-42720 + J (7 + J) (-104 + J (7 + J))) mA^10 s^2 - 640 (-3 + J) (-1 + J) (8 + J) (10 + J) (-12480 + J (7 + J) (-104 + J (7 + J))) mA^8 s^3 + 320 (-3 + J) (-1 + J) (8 + J) (10 + J) (-4920 + J (7 + J) (-104 + J (7 + J))) mA^6 s^4 - 80 (-3 + J) (-1 + J) (8 + J) (10 + J) (-1896 + J (7 + J) (-104 + J (7 + J))) mA^4 s^5 + 10 (-3 + J) (-1 + J) (8 + J) (10 + J) (-48 + (-1 + J) J) (8 + J (15 + J)) mA^2 s^6 - (35064000 + J (7 + J) (-1763640 + J (7 + J) (31942 + J (7 + J) (-277 + J (7 + J))))) s^7);
  N[result, 600]
];

x110snd = Module[{s, mA, J, result},
  s = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = 1/(-4 mA^2 + s)^12 - 1/(289989867847680000 (-4 mA^2 + s)^12) (-10 + J) (-8 + J) (-6 + J) (-4 + J) (-2 + J) J (7 + J) (9 + J) (11 + J) (13 + J) (15 + J) (17 + J) (-833676480 + J (7 + J) (114378912 + J (7 + J) (-4230052 + J (7 + J) (63448 + J (7 + J) (-417 + J (7 + J)))))) + (16777216 mA^24)/(s^12 (-4 mA^2 + s)^12) - (50331648 mA^22)/(s^11 (-4 mA^2 + s)^12) + (69206016 mA^20)/(s^10 (-4 mA^2 + s)^12) - (57671680 mA^18)/(s^9 (-4 mA^2 + s)^12) + (32440320 mA^16)/(s^8 (-4 mA^2 + s)^12) - (12976128 mA^14)/(s^7 (-4 mA^2 + s)^12) + (3784704 mA^12)/(s^6 (-4 mA^2 + s)^12) - (811008 mA^10)/(s^5 (-4 mA^2 + s)^12) + (126720 mA^8)/(s^4 (-4 mA^2 + s)^12) - (14080 mA^6)/(s^3 (-4 mA^2 + s)^12) + (1056 mA^4)/(s^2 (-4 mA^2 + s)^12) + 1/(144994933923840000 s (-4 mA^2 + s)^12) (-6959756828344320000 + (-10 + J) (-9 + J) (-8 + J) (-7 + J) (-6 + J) (-5 + J) (-4 + J) (-3 + J) (-2 + J) (-1 + J) J (7 + J) (8 + J) (9 + J) (10 + J) (11 + J) (12 + J) (13 + J) (14 + J) (15 + J) (16 + J) (17 + J)) mA^2;
  N[result, 600]
]; *)

(* need simplification *)
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
