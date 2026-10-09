(* ::Package:: *)

sourceDirectory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{sourceDirectory, "..", "SDPB.m"}]];

(* Keep the spectral data exact until the final numerical conversion. *)
m1 = 55/100;
J1 = 0;
J2 = 2;
mgap = 83/50;

nulllist = {30, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];
functionalDimension = 2 + Length[list0];
massPrefactorPower = 18;
JCasimirPower = 9;

Nlist[n_, z_, J_] := {z-1/2 J (7+J) z,1-1/20 J (7+J) (-13+J (7+J)),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z)),1/z^2-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^2),-((J (7+J) (-23+J (7+J)))/(20 z^2)),1/z^3-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)))/(403200 z^3),-((J (7+J) (540+J (7+J) (-53+J (7+J))))/(360 z^3)),1/z^4-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J)))/(21772800 z^4),-((J (7+J) (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))))/(10080 z^4)),1/z^5-((-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J))))/(1524096000 z^5),-((J (7+J) (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))))/(403200 z^5)),-(((-2+J) J (7+J) (9+J) (-53+J (7+J)))/(360 z^5)),1/z^6-((-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J))))))/(134120448000 z^6),-((J (7+J) (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))))/(21772800 z^6)),-(((-2+J) J (7+J) (9+J) (2564+J (7+J) (-108+J (7+J))))/(5040 z^6)),1/z^7-((-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J))))))/(14485008384000 z^7),-((J (7+J) (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))))/(1524096000 z^7)),-(((-2+J) J (7+J) (9+J) (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))))/(201600 z^7)),1/z^8-1/(1883051089920000 z^8)(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J)))))),-1/(134120448000 z^8)J (7+J) (-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))),-(((-2+J) J (7+J) (9+J) (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))))/(10886400 z^8)),-(((-2+J) J (7+J) (9+J) (35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))))/(21772800 z^8)),1/z^9-1/(289989867847680000 z^9)(-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J)))))),-1/(14485008384000 z^9)J (7+J) (180894269952000+J (7+J) (-33122195274240+J (7+J) (2349433112832+J (7+J) (-85849259040+J (7+J) (1788748560+J (7+J) (-22054832+J (7+J) (158952+J (7+J) (-618+J (7+J))))))))),-(((-2+J) J (7+J) (9+J) (-2563473600+J (7+J) (175893840+J (7+J) (-4726576+J (7+J) (61658+J (7+J) (-395+J (7+J)))))))/(762048000 z^9)),-(((-2+J) J (7+J) (9+J) (-5573865600+J (7+J) (318324240+J (7+J) (-6824476+J (7+J) (71108+J (7+J) (-395+J (7+J)))))))/(1524096000 z^9)),1/z^10-1/(52198176212582400000 z^10)(-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (172008748800+J (7+J) (-25008497280+J (7+J) (1017235728+J (7+J) (-17785440+J (7+J) (152128+J (7+J) (-628+J (7+J))))))),-1/(1883051089920000 z^10)J (7+J) (-30370283513856000+J (7+J) (5685348310272000+J (7+J) (-431278994177280+J (7+J) (17117281051776+J (7+J) (-396920087856+J (7+J) (5654064144+J (7+J) (-50058872+J (7+J) (268176+J (7+J) (-795+J (7+J)))))))))),-1/(67060224000 z^10)(-2+J) J (7+J) (9+J) (357678604800+J (7+J) (-27684313920+J (7+J) (855078608+J (7+J) (-13637120+J (7+J) (118668+J (7+J) (-538+J (7+J))))))),-1/(134120448000 z^10)(-2+J) J (7+J) (9+J) (942637132800+J (7+J) (-59261829120+J (7+J) (1465072608+J (7+J) (-18734520+J (7+J) (134068+J (7+J) (-538+J (7+J))))))),1/z^11-1/(10857220652217139200000 z^11)(-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (178803072000+J (7+J) (-26522288640+J (7+J) (1083860832+J (7+J) (-18822392+J (7+J) (158628+J (7+J) (-642+J (7+J)))))))}[[n+1]];







(* ---------------------------------------------------------------------- *)
(* Functional vectors and polynomial blocks                               *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  rawFunctionalVector,
  polyify,
  polynomializeVector,
  JScale,
  makePolynomialBlock,
  continuumBlock,
  fixedStateBlock
];

rawFunctionalVector[z_, J_, secondEntry_: 0] := Join[
  {2z^2, secondEntry},
  Table[Nlist[n, z, J], {n, 0, nulllist[[1]]}]
];

polyify[expr_, var_Symbol] := Expand[Cancel[Together[expr]]];

polynomializeVector[vector_List, z_Symbol] :=
  polyify[z^massPrefactorPower #, z] & /@ vector;

(* A positive blockwise rescaling leaves the feasible cone unchanged, but
   removes the leading q^9 = O(J^18) growth, where q = J (J + 7). *)
JScale[J_Integer] := 1/(1 + J (J + 7))^JCasimirPower;

makePolynomialBlock[
  J_Integer,
  zValue_,
  variable_Symbol,
  secondEntry_: 0,
  useJScaling_: True
] := Module[{zInternal, vector, badComponents, scale},
  vector = polynomializeVector[
    rawFunctionalVector[zInternal, J, secondEntry],
    zInternal
  ];

  badComponents = Flatten @ Position[
    PolynomialQ[#, zInternal] & /@ vector,
    False
  ];
  If[badComponents =!= {},
    Print["Non-polynomial components before the mass shift: ", badComponents];
    Abort[]
  ];

  vector = Expand[# /. zInternal -> zValue] & /@ vector;
  badComponents = Flatten @ Position[
    PolynomialQ[#, variable] & /@ vector,
    False
  ];
  If[badComponents =!= {},
    Print["Non-polynomial components after the mass shift: ", badComponents];
    Abort[]
  ];

  If[Length[vector] =!= functionalDimension,
    Print[
      "Functional dimension mismatch: expected ", functionalDimension,
      ", received ", Length[vector]
    ];
    Abort[]
  ];

  scale = If[TrueQ[useJScaling], JScale[J], 1];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    {{scale vector}}
  ]
];

continuumBlock[J_Integer, variable_Symbol] :=
  makePolynomialBlock[J, mgap + variable, variable, 0, True];

fixedStateBlock[
  J_Integer,
  massSquared_,
  variable_Symbol,
  secondEntry_: 0
] := makePolynomialBlock[
  J,
  massSquared,
  variable,
  secondEntry,
  False
];

(* ---------------------------------------------------------------------- *)
(* Analytically generated asymptotic constraints                          *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  fixedMassLargeJVector,
  fixedMassLargeJBlock,
  toCasimir,
  impactVector,
  impactBlock
];

(* Fixed mass, J -> Infinity.  Dividing by the largest power of J is a
   positive rescaling and reproduces the old hand-written Polyinf block. *)
fixedMassLargeJVector[z_] := Module[
  {J, vector, rationalVector, degrees, maxDegree},

  vector = rawFunctionalVector[z, J, 0];
  rationalVector = Together /@ vector;
  degrees = Exponent[Numerator[#], J] & /@ rationalVector;
  maxDegree = Max[degrees];

  MapThread[
    Function[{entry, degree},
      If[
        degree === maxDegree,
        Coefficient[Numerator[entry], J, maxDegree]/Denominator[entry],
        0
      ]
    ],
    {rationalVector, degrees}
  ]
];

fixedMassLargeJBlock[variable_Symbol] := Module[
  {zInternal, vector, badComponents},

  vector = polynomializeVector[
    fixedMassLargeJVector[zInternal],
    zInternal
  ];
  vector = Expand[# /. zInternal -> mgap + variable] & /@ vector;

  badComponents = Flatten @ Position[
    PolynomialQ[#, variable] & /@ vector,
    False
  ];
  If[badComponents =!= {},
    Print["Invalid fixed-mass large-J components: ", badComponents];
    Abort[]
  ];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    {{vector}}
  ]
];

(* Every Nlist entry is invariant under J -> -7-J and can therefore be
   written in terms of the ten-dimensional Casimir q = J(J+7). *)
toCasimir[expr_, J_Symbol, casimir_Symbol] := Module[
  {entry, numerator, denominator, remainder},

  entry = Together[expr];
  numerator = Numerator[entry];
  denominator = Denominator[entry];

  remainder = PolynomialRemainder[
    numerator,
    J^2 + 7 J - casimir,
    J
  ];

  If[!FreeQ[remainder, J] || !FreeQ[denominator, J],
    Print["Could not rewrite a null constraint in terms of J(J+7)."];
    Abort[]
  ];

  Cancel[remainder/denominator]
];

(* Correlated large-mass/large-J limit with r = J(J+7)/z fixed.
   This is the missing impact-parameter-type boundary of the moment cone. *)
impactVector[r_] := Module[
  {z, J, casimir, vector, casimirVector, result},

  vector = rawFunctionalVector[z, J, 0];
  casimirVector = toCasimir[#, J, casimir] & /@ vector;

  result = FullSimplify[
    Limit[
      (casimirVector /. casimir -> r z)/z^2,
      z -> Infinity
    ],
    Assumptions -> r >= 0
  ];

  If[
    !FreeQ[result, DirectedInfinity | Indeterminate] ||
    !AllTrue[result, PolynomialQ[#, r] &],
    Print["The correlated large-J limit is not a finite polynomial vector."];
    Abort[]
  ];

  result
];

impactBlock[variable_Symbol] := Module[{vector},
  vector = impactVector[variable];

  If[Length[vector] =!= functionalDimension,
    Print["Impact-vector dimension mismatch."];
    Abort[]
  ];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    {{vector}}
  ]
];

(* ---------------------------------------------------------------------- *)
(* PMP assembly                                                           *)
(* ---------------------------------------------------------------------- *)

(* Dense low J values locate the usual extremal states.  The sparse probes
   monitor the transition to the analytic asymptotic blocks.  A final
   solution must still be checked for omitted finite-J violations. *)
coreJMax = 200;
JProbeList = {250, 300, 400, 500, 700, 1000, 1500, 2500, 5000};
continuumJList = DeleteDuplicates @ Join[
  Range[0, coreJMax, 2],
  JProbeList
];

(* Post-SDPB certification helper.  Supply the dual functional in the same
   component order as obj/norm.  Any returned J must be added to
   continuumJList before regenerating the PMP. *)
ClearAll[
  continuumPolynomialVector,
  minimumAtJ,
  findViolatingJs
];

continuumPolynomialVector[J_Integer, variable_Symbol] := Module[
  {zInternal, vector},

  vector = polynomializeVector[
    rawFunctionalVector[zInternal, J, 0],
    zInternal
  ];

  JScale[J] (Expand[# /. zInternal -> mgap + variable] & /@ vector)
];

minimumAtJ[
  functional_List,
  J_Integer,
  variable_Symbol,
  prec_: 80
] := Module[{polynomial},
  If[Length[functional] =!= functionalDimension,
    Print["Dual-functional dimension mismatch."];
    Abort[]
  ];

  polynomial = N[
    functional . continuumPolynomialVector[J, variable],
    prec
  ];

  NMinimize[
    {polynomial, variable >= 0},
    variable,
    WorkingPrecision -> prec,
    AccuracyGoal -> Floor[prec/3],
    PrecisionGoal -> Floor[prec/3]
  ]
];

findViolatingJs[
  functional_List,
  JValues_List,
  tolerance_: 10^-30,
  prec_: 80
] := DeleteCases[
  Table[
    Module[{minimum},
      minimum = Quiet @ Check[
        minimumAtJ[functional, J, x, prec],
        $Failed
      ];

      Which[
        minimum === $Failed,
          <|"J" -> J, "Status" -> "MinimizationFailed"|>,
        First[minimum] < -Abs[tolerance],
          <|
            "J" -> J,
            "Minimum" -> First[minimum],
            "Location" -> (x /. Last[minimum])
          |>,
        True,
          Nothing
      ]
    ],
    {J, JValues}
  ],
  Nothing
];

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
  {
    continuumBlocks, specialBlocks, asymptoticBlocks,
    pols, norm, obj, expectedBlocks
  },

  If[!TrueQ[0 < m1 < 1 < mgap],
    Print["Invalid spectrum ordering: expected 0 < m1 < 1 < mgap."];
    Abort[]
  ];

  Print["Building ", Length[continuumJList], " finite-J blocks..."];
  DistributeDefinitions[
    Nlist,
    nulllist,
    massPrefactorPower,
    JCasimirPower,
    functionalDimension,
    mgap,
    rawFunctionalVector,
    polyify,
    polynomializeVector,
    JScale,
    makePolynomialBlock,
    continuumBlock
  ];
  continuumBlocks = ParallelMap[
    continuumBlock[#, x] &,
    continuumJList
  ];

  specialBlocks = {
    fixedStateBlock[J2, 1, x, 1],
    fixedStateBlock[J1, m1, x, 0]
  };

  Print["Building analytic large-J blocks..."];
  asymptoticBlocks = {
    fixedMassLargeJBlock[x],
    impactBlock[x]
  };

  pols = N[
    Join[specialBlocks, continuumBlocks, asymptoticBlocks],
    prec
  ];

  expectedBlocks = 4 + Length[continuumJList];
  If[Length[pols] =!= expectedBlocks,
    Print[
      "Unexpected block count: built ", Length[pols],
      ", expected ", expectedBlocks, "."
    ];
    Abort[]
  ];

  norm = -N[Flatten[{{0, 1}, list0}], prec];
  obj = -N[Flatten[{{1, 0}, list0}], prec];

  If[
    Length[norm] =!= functionalDimension ||
    Length[obj] =!= functionalDimension,
    Print["Objective or normalization dimension mismatch."];
    Abort[]
  ];

  Print["functional dimension = ", functionalDimension];
  Print["number of PMP blocks = ", Length[pols]];
  Print["finite J values = ", continuumJList];

  WritePmpJson[
    datfile,
    SDP[obj, norm, pols],
    prec,
    getAnalyticSampleData
  ]
];

outputFile = FileNameJoin[{sourceDirectory, "n_pmp.json"}];
If[!TrueQ[$Test16SkipExport],
  PMP2SDP[outputFile, 1000]
];
