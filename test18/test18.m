sourceDirectory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{sourceDirectory, "..", "SDPB.m"}]];

(* Keep physical input exact until the numerical blocks are assembled. *)
m1 = 2/5;
mA = 1/1000;
J1 = 0;
J2 = 2;
mgap = 83/50;

(* 7 null constraints *)
nulllist = {6, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];

NBBBB[n_, z_, J_] := {z-1/2 J (7+J) z,1-1/20 J (7+J) (-13+J (7+J)),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z)),1/z^2-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^2),-((J (7+J) (-23+J (7+J)))/(20 z^2)),1/z^3-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)))/(403200 z^3),-((J (7+J) (540+J (7+J) (-53+J (7+J))))/(360 z^3))}[[n+1]];


NAAAA[n_, x_, J_]:= {((-4 mA^2+x)^(3/2) (32 mA^4+2 (-1+J) (8+J) mA^2 x-(-2+J (7+J)) x^2))/(2 x^(5/2)),1/(20 x^(7/2))Sqrt[-4 mA^2+x] (-1280 mA^6+960 mA^4 x+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 x^2-(-20+J (7+J) (-13+J (7+J))) x^3),1/(360 x^(9/2) Sqrt[-4 mA^2+x])(92160 mA^8-92160 mA^6 x+34560 mA^4 x^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 x^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) x^4),-((J (7+J) mA^2 (80 mA^4+2 (-38+J (7+J)) mA^2 x-(-23+J (7+J)) x^2))/(5 x^(7/2) Sqrt[-4 mA^2+x])),1/(10080 x^(11/2) (-4 mA^2+x)^(3/2))(-10321920 mA^10+12902400 mA^8 x-6451200 mA^6 x^2+1612800 mA^4 x^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 x^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) x^5),1/(180 x^(9/2) (-4 mA^2+x)^(3/2))J (7+J) (11520 mA^8-11520 mA^6 x-4 (-936+J (7+J) (-26+J (7+J))) mA^4 x^2+2 (-216+J (7+J) (-26+J (7+J))) mA^2 x^3-9 (-23+J (7+J)) x^4),1/(403200 x^(13/2) (-4 mA^2+x)^(5/2))(1651507200 mA^12-2477260800 mA^10 x+1548288000 mA^8 x^2-516096000 mA^6 x^3+96768000 mA^4 x^4+2 (-4838400+(-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J)) mA^2 x^5-(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) x^6)}[[n+1]];


(* ---------------------------------------------------------------------- *)
(* Analytic half-line PMP formulation, matching test16.m                  *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  phaseAAAA,
  phaseBBBB,
  rawAAAA,
  rawBBBB,
  rawChannelVector,
  validatePhysicalDomain,
  channelClearingFactor,
  polynomializeChannelVector,
  spinScale,
  makeChannelBlock,
  continuumBlocksForSpin,
  fixedStateBlocks,
  jDegree,
  fixedMassLargeJVector,
  fixedMassLargeJBlock,
  toCasimir,
  impactVector,
  lowestPolynomialPower,
  stripCommonBoundaryPower,
  impactBlock
];

functionalCount = 2 + Length[list0];
channelNames = {"AAAA", "BBBB"};

validatePhysicalDomain[] := If[
  !TrueQ[mgap > 4 mA^2 && m1 > 4 mA^2 && 1 > 4 mA^2],
  Print["All continuum and isolated-state masses must satisfy z > 4 mA^2."];
  Abort[]
];

(* The following factors are strictly positive on every continuum and
   isolated-state point used below.  Dividing a diagonal channel by one of
   them is therefore an invertible congruence rescaling and preserves the
   positivity problem exactly. *)
phaseAAAA[z_] := (z - 4 mA^2)^(7/2)/Sqrt[z];
phaseBBBB[z_] := z^3;

(* PowerExpand is safe here because every use satisfies z > 4 mA^2 > 0. *)
rawAAAA[n_Integer, z_, spin_] := Cancel @ Together @ PowerExpand[
  NAAAA[n, z, spin]/phaseAAAA[z]
];

rawBBBB[n_Integer, z_, spin_] := Cancel @ Together[
  NBBBB[n, z, spin]/phaseBBBB[z]
];

(* Coefficient order: g0, normalized J=2 coupling, seven null directions. *)
rawChannelVector[
  "AAAA",
  z_,
  spin_,
  normalizedEntry_: 0
] := Join[
  {0, normalizedEntry},
  Table[rawAAAA[n, z, spin], {n, 0, nulllist[[1]]}]
];

rawChannelVector[
  "BBBB",
  z_,
  spin_,
  normalizedEntry_: 0
] := Join[
  {2/z, normalizedEntry},
  Table[rawBBBB[n, z, spin], {n, 0, nulllist[[1]]}]
];

(* The first seven kernels require at most six powers of each physical
   denominator.  These factors are positive throughout the relevant domain. *)
channelClearingFactor["AAAA", z_] := z^6 (z - 4 mA^2)^6;
channelClearingFactor["BBBB", z_] := z^6;

polynomializeChannelVector[channel_String, vector_List, z_Symbol] :=
  Expand[Cancel[Together[channelClearingFactor[channel, z] #]]] & /@ vector;

(* The largest power in spin is J^10 = O[(J(J+7))^5]. *)
spinScalePower = 5;
spinScale[spin_Integer] := 1/(1 + spin (spin + 7))^spinScalePower;

makeChannelBlock[
  channel_String,
  spin_Integer,
  zValue_,
  variable_Symbol,
  normalizedEntry_: 0,
  useSpinScaling_: True
] := Module[{zInternal, vector, badComponents, scale},
  vector = polynomializeChannelVector[
    channel,
    rawChannelVector[channel, zInternal, spin, normalizedEntry],
    zInternal
  ];
  vector = Expand[# /. zInternal -> zValue] & /@ vector;

  badComponents = Flatten @ Position[
    PolynomialQ[#, variable] & /@ vector,
    False
  ];
  If[badComponents =!= {},
    Print[
      "Non-polynomial ", channel, " components at spin ", spin,
      ": ", badComponents
    ];
    Abort[]
  ];

  If[Length[vector] =!= functionalCount,
    Print[
      "Functional dimension mismatch in ", channel, ": expected ",
      functionalCount, ", received ", Length[vector], "."
    ];
    Abort[]
  ];

  scale = If[TrueQ[useSpinScaling], spinScale[spin], 1];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    {{scale vector}}
  ]
];

continuumBlocksForSpin[spin_Integer, variable_Symbol] :=
  makeChannelBlock[#, spin, mgap + variable, variable, 0, True] & /@
    channelNames;

fixedStateBlocks[
  spin_Integer,
  massSquared_,
  variable_Symbol,
  normalizedEntry_: 0
] := makeChannelBlock[
  #,
  spin,
  massSquared,
  variable,
  normalizedEntry,
  True
] & /@ channelNames;

(* ---------------------------------------------------------------------- *)
(* Analytic high-spin boundaries                                          *)
(* ---------------------------------------------------------------------- *)

jDegree[expr_, spin_Symbol] := Module[{rational},
  If[TrueQ[expr === 0],
    -Infinity,
    rational = Together[expr];
    Exponent[Numerator[rational], spin] -
      Exponent[Denominator[rational], spin]
  ]
];

fixedMassLargeJVector[channel_String, z_] := Module[
  {spin, vector, rationalVector, degrees, maxDegree},
  vector = rawChannelVector[channel, z, spin, 0];
  rationalVector = Together /@ vector;
  degrees = jDegree[#, spin] & /@ rationalVector;
  maxDegree = Max[degrees];

  MapThread[
    Function[{entry, degree},
      If[
        degree === maxDegree,
        Coefficient[Numerator[entry], spin, maxDegree]/Denominator[entry],
        0
      ]
    ],
    {rationalVector, degrees}
  ]
];

fixedMassLargeJBlock[channel_String, variable_Symbol] := Module[
  {zInternal, vector, badComponents},
  vector = polynomializeChannelVector[
    channel,
    fixedMassLargeJVector[channel, zInternal],
    zInternal
  ];
  vector = Expand[# /. zInternal -> mgap + variable] & /@ vector;

  badComponents = Flatten @ Position[
    PolynomialQ[#, variable] & /@ vector,
    False
  ];
  If[badComponents =!= {},
    Print["Invalid fixed-mass large-J block for ", channel, "."];
    Abort[]
  ];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    {{vector}}
  ]
];

toCasimir[expr_, spin_Symbol, casimir_Symbol] := Module[
  {entry, numerator, denominator, remainder},
  entry = Together[expr];
  numerator = Numerator[entry];
  denominator = Denominator[entry];
  remainder = PolynomialRemainder[
    numerator,
    spin^2 + 7 spin - casimir,
    spin
  ];

  If[!FreeQ[remainder, spin] || !FreeQ[denominator, spin],
    Print["Could not rewrite a kernel in J(J+7)."];
    Abort[]
  ];

  Cancel[remainder/denominator]
];

(* Correlated boundary z,J -> Infinity with r=J(J+7)/z fixed. *)
impactVector[channel_String, r_] := Module[
  {z, spin, casimir, vector, casimirVector, result},
  vector = rawChannelVector[channel, z, spin, 0];
  casimirVector = toCasimir[#, spin, casimir] & /@ vector;
  result = FullSimplify[
    Limit[z (casimirVector /. casimir -> r z), z -> Infinity],
    Assumptions -> r >= 0
  ];

  If[
    !FreeQ[result, DirectedInfinity[___] | ComplexInfinity | Indeterminate] ||
    !AllTrue[result, PolynomialQ[#, r] &],
    Print["Invalid correlated large-spin vector for ", channel, "."];
    Abort[]
  ];

  result
];

(* If every component of an impact vector contains the same power r^k,
   remove it before constructing the PMP block.  For r >= 0 and polynomial
   q(r), r^k q(r) >= 0 is equivalent to q(r) >= 0: the statement is immediate
   for r > 0, while the endpoint follows by continuity.  Keeping a common
   zero creates an artificial null interpolation direction and a singular
   Schur block in SDPB. *)
lowestPolynomialPower[expr_, variable_Symbol] := Module[{rules},
  If[TrueQ[PossibleZeroQ[expr]], Return[Infinity]];
  rules = CoefficientRules[Expand[expr], variable];
  If[rules === {}, Infinity, Min[First /@ rules[[All, 1]]]]
];

stripCommonBoundaryPower[vector_List, variable_Symbol] := Module[
  {orders, commonPower, reducedVector},
  orders = DeleteCases[
    lowestPolynomialPower[#, variable] & /@ vector,
    Infinity
  ];
  commonPower = If[orders === {}, 0, Min[orders]];
  reducedVector = Expand[
    Cancel[Together[#/variable^commonPower]]
  ] & /@ vector;
  {reducedVector, commonPower}
];

impactBlock[channel_String, variable_Symbol] := Module[
  {rawVector, vector, commonPower},
  rawVector = impactVector[channel, variable];
  {vector, commonPower} = stripCommonBoundaryPower[rawVector, variable];

  If[commonPower > 0,
    Print[
      "Removed the common boundary factor ", variable, "^", commonPower,
      " from the ", channel, " impact block."
    ]
  ];

  If[Length[vector] =!= functionalCount,
    Print["Impact-vector dimension mismatch for ", channel, "."];
    Abort[]
  ];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    {{vector}}
  ]
];

(* Match the spin treatment of test16.m. *)
coreSpinMax = 200;
spinProbeList = {250, 300, 400, 500, 700, 1000, 1500, 2500, 5000};
continuumSpinList = DeleteDuplicates @ Join[
  Range[0, coreSpinMax, 2],
  spinProbeList
];

(* Keep the zero-objective feasibility test as the default.  Once this run
   converges, set feasibilityOnly=False and choose objectiveSign=+1 or -1. *)
feasibilityOnly = True;
objectiveSign = -1;

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
  {continuumBlocks, specialBlocks, asymptoticBlocks, pols, norm, obj},

  validatePhysicalDomain[];

  DistributeDefinitions[
    mA, mgap, nulllist, list0, functionalCount, channelNames,
    NBBBB, NAAAA, phaseAAAA, phaseBBBB, rawAAAA, rawBBBB,
    rawChannelVector, channelClearingFactor, polynomializeChannelVector,
    spinScalePower, spinScale, makeChannelBlock, continuumBlocksForSpin
  ];

  Print[
    "Building analytic continuum blocks for ",
    Length[continuumSpinList], " spins and two channels..."
  ];
  continuumBlocks = Flatten[
    ParallelMap[continuumBlocksForSpin[#, x] &, continuumSpinList],
    1
  ];

  specialBlocks = Join[
    fixedStateBlocks[J1, m1, x, 0],
    fixedStateBlocks[J2, 1, x, 1]
  ];

  Print["Building fixed-mass and correlated asymptotic blocks..."];
  asymptoticBlocks = Join[
    fixedMassLargeJBlock[#, x] & /@ channelNames,
    impactBlock[#, x] & /@ channelNames
  ];

  pols = N[
    Join[specialBlocks, continuumBlocks, asymptoticBlocks],
    prec
  ];

  If[Length[pols] =!= 2 Length[continuumSpinList] + 8,
    Print["Unexpected number of PMP blocks: ", Length[pols], "."];
    Abort[]
  ];

  norm = -N[Flatten[{{0, 1}, list0}], prec];
  obj = If[
    TrueQ[feasibilityOnly],
    ConstantArray[0, functionalCount],
    objectiveSign N[Flatten[{{1, 0}, list0}], prec]
  ];

  If[Length[norm] =!= functionalCount || Length[obj] =!= functionalCount,
    Print["Objective or normalization dimension mismatch."];
    Abort[]
  ];

  Print["functional dimension = ", functionalCount];
  Print["finite spins = ", continuumSpinList];
  Print["PMP blocks = ", Length[pols]];
  Print["feasibility-only mode = ", feasibilityOnly];
  Print["Writing ", datfile, "..."];

  WritePmpJson[
    datfile,
    SDP[obj, norm, pols],
    prec,
    getAnalyticSampleData
  ];

  Print["Wrote ", datfile, "."]
];

outputFile = FileNameJoin[{sourceDirectory, "n_pmp.json"}];
If[!TrueQ[$Test18SkipExport],
  PMP2SDP[outputFile, 600]
];
