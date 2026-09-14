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


Nlist[n_,z_,J_] := {
	{NAAAA[n,z,J],0},
	{0,NBBBB[n,z,J]}
};


polyify[expr_] := Expand @ Cancel @ Together[expr];

PolyInfBBBB[n_,J_,x_] := {0,0,0,0,0,-(1/(403200 x^3)),0}[[n+1]];


PolyInfAAAA[n_,J_,x_] := {0,0,0,0,0,0,(2 mA^2-x)/(403200 x^(3/2) (-4 mA^2+x)^(5/2))}[[n+1]];


NPolyInf[n_,J_,x_] := {
  {PolyInfAAAA[n,J,x],0},
  {0,PolyInfBBBB[n,J,x]}
};

(* ---------------------------------------------------------------------- *)
(* Paper-inspired sampling and numerical conditioning                     *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  chebyshevPhi,
  compactCoordinate,
  energyFromPhi,
  paperSamples,
  zeroMatrix,
  g0CoordinateMatrix,
  lambda22CoordinateMatrix,
  coefficientMatrices,
  largeJCoefficientMatrices,
  congruenceRescale,
  coefficientTensor,
  pmpBlock,
  finiteStateBlock,
  largeJBlock,
  validateConfiguration,
  validateSampling
];

matrixDimension = 2;
functionalCount = 2 + Length[list0];

(* Appendix D, eq. (169): Chebyshev nodes in the conformal angle phi.
   The endpoints phi=0 and phi=Pi are intentionally not sampled. *)
nPoints = 200;
lMax = 60;
zReference = 0;
useCongruenceRescaling = True;

chebyshevPhi[k_Integer, count_Integer] :=
  Pi/2 + Pi/2 Cos[((k + 1/2) Pi)/count];

(* For rho(z,mgap,zReference)=Exp[I phi], inversion of the conformal map gives
     z = mgap + (mgap-zReference) Tan[phi/2]^2.
   With zReference=0 this is equivalent to the old map
     z=mgap/(1-x), x=Sin[phi/2]^2. *)
compactCoordinate[phi_] := Sin[phi/2]^2;

energyFromPhi[phi_] :=
  mgap + (mgap - zReference) Tan[phi/2]^2;

paperSamples[count_Integer, prec_Integer] := SortBy[
  Table[
    With[{phi = N[chebyshevPhi[k, count], prec]},
      <|
        "Index" -> k,
        "Phi" -> phi,
        "CompactCoordinate" -> N[compactCoordinate[phi], prec],
        "Energy" -> N[energyFromPhi[phi], prec]
      |>
    ],
    {k, 0, count - 1}
  ],
  #["Energy"] &
];

zeroMatrix = ConstantArray[0, {matrixDimension, matrixDimension}];

(* First functional coordinate: g0. *)
g0CoordinateMatrix[z_] := {
  {0, 0},
  {0, 2 z^2}
};

(* Second functional coordinate: the on-shell coupling of the fixed J=2 state. *)
lambda22CoordinateMatrix[z_] := {
  {(-4 mA^2 + z)^(7/2)/Sqrt[z], 0},
  {0, z^3}
};

coefficientMatrices[
  z_,
  spin_Integer,
  lambda22Matrix_: zeroMatrix
] := Join[
  {g0CoordinateMatrix[z], lambda22Matrix},
  Table[Nlist[n, z, spin], {n, 0, nulllist[[1]]}]
];

largeJCoefficientMatrices[z_] := Join[
  {zeroMatrix, zeroMatrix},
  Table[NPolyInf[n, 0, z], {n, 0, nulllist[[1]]}]
];

(* Paper eqs. (174)-(175), adapted to the two diagonal channels used here.
   Every coefficient matrix of one constraint receives the same congruence
   transformation, so the PSD cone and the bound are unchanged. *)
congruenceRescale[matrices_List, prec_Integer] := Module[
  {numericMatrices, channelNorms, diagonalRescaling},

  numericMatrices = N[matrices, prec];
  channelNorms = Table[
    Max @@ Abs[numericMatrices[[All, channel, channel]]],
    {channel, matrixDimension}
  ];
  channelNorms = Replace[channelNorms, value_ /; TrueQ[value == 0] -> 1, {1}];

  diagonalRescaling = DiagonalMatrix[1/Sqrt[channelNorms]];
  N[diagonalRescaling . # . diagonalRescaling, prec] & /@ numericMatrices
];

coefficientTensor[matrices_List] := Table[
  Table[
    matrices[[coefficient, row, column]],
    {coefficient, Length[matrices]}
  ],
  {row, matrixDimension},
  {column, matrixDimension}
];

pmpBlock[matrices_List, variable_Symbol, prec_Integer] := Module[
  {preparedMatrices, tensor},

  If[Dimensions[matrices] =!= {functionalCount, matrixDimension, matrixDimension},
    Print[
      "Invalid coefficient-matrix dimensions: ", Dimensions[matrices],
      "; expected ", {functionalCount, matrixDimension, matrixDimension}, "."
    ];
    Abort[]
  ];

  If[matrices =!= Transpose[matrices, {1, 3, 2}],
    Print["A sampled coefficient matrix is not symmetric."];
    Abort[]
  ];

  preparedMatrices = If[
    TrueQ[useCongruenceRescaling],
    congruenceRescale[matrices, prec],
    N[matrices, prec]
  ];
  tensor = coefficientTensor[preparedMatrices];

  If[!FreeQ[tensor, Indeterminate | ComplexInfinity | DirectedInfinity[___]],
    Print["A sampled block contains a non-finite coefficient."];
    Abort[]
  ];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, variable],
    tensor
  ]
];

finiteStateBlock[
  z_,
  spin_Integer,
  variable_Symbol,
  lambda22Matrix_: zeroMatrix,
  prec_: 600
] := pmpBlock[
  coefficientMatrices[z, spin, lambda22Matrix],
  variable,
  prec
];

largeJBlock[z_, variable_Symbol, prec_: 600] := pmpBlock[
  largeJCoefficientMatrices[z],
  variable,
  prec
];

validateConfiguration[] := Module[{},
  If[!IntegerQ[nPoints] || !(1 <= nPoints <= 200),
    Print["nPoints must be an integer between 1 and 200."];
    Abort[]
  ];

  If[!IntegerQ[lMax] || lMax < 0 || OddQ[lMax],
    Print["lMax must be a non-negative even integer."];
    Abort[]
  ];

  If[!TrueQ[zReference < mgap],
    Print["zReference must lie below the continuum threshold mgap."];
    Abort[]
  ];
];

validateSampling[samples_List, spins_List] := Module[
  {phis, energies, compactValues},

  phis = Lookup[samples, "Phi"];
  energies = Lookup[samples, "Energy"];
  compactValues = Lookup[samples, "CompactCoordinate"];

  If[
    Length[samples] =!= nPoints ||
    Length[DeleteDuplicates[energies]] =!= nPoints,
    Print["The energy sample contains missing or duplicate points."];
    Abort[]
  ];

  If[
    !AllTrue[phis, 0 < # < Pi &] ||
    !AllTrue[compactValues, 0 < # < 1 &] ||
    !AllTrue[energies, # > mgap &],
    Print["A Chebyshev sample lies outside the physical domain."];
    Abort[]
  ];

  If[spins =!= Range[0, lMax, 2],
    Print["Spin list is inconsistent with lMax."];
    Abort[]
  ];
];

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
  {
    samples, spins, continuumBlocks, largeJBlocks,
    specialBlocks, pols, norm, obj
  },

  validateConfiguration[];
  samples = paperSamples[nPoints, prec];
  spins = Range[0, lMax, 2];
  validateSampling[samples, spins];

  Print["Chebyshev energy samples = ", Length[samples]];
  Print["phi range = ", {First[samples]["Phi"], Last[samples]["Phi"]}];
  Print["energy range = ", {First[samples]["Energy"], Last[samples]["Energy"]}];
  Print["even spins = 0, 2, ..., ", lMax, " (", Length[spins], " spins)"];

  DistributeDefinitions[
    mA, mgap, nulllist, list0, matrixDimension, functionalCount,
    NBBBB, NAAAA, Nlist, PolyInfBBBB, PolyInfAAAA, NPolyInf,
    zeroMatrix, g0CoordinateMatrix, lambda22CoordinateMatrix,
    coefficientMatrices, largeJCoefficientMatrices,
    useCongruenceRescaling, congruenceRescale, coefficientTensor,
    pmpBlock, finiteStateBlock, largeJBlock
  ];

  Print["Building finite-spin sampled blocks..."];
  continuumBlocks = Flatten[
    ParallelTable[
      finiteStateBlock[sample["Energy"], spin, x, zeroMatrix, prec],
      {sample, samples},
      {spin, spins}
    ],
    1
  ];

  Print["Building one analytic large-J block per energy sample..."];
  largeJBlocks = ParallelMap[
    largeJBlock[#["Energy"], x, prec] &,
    samples
  ];

  specialBlocks = {
    finiteStateBlock[m1, J1, x, zeroMatrix, prec],
    finiteStateBlock[1, J2, x, lambda22CoordinateMatrix[1], prec]
  };

  pols = Join[specialBlocks, continuumBlocks, largeJBlocks];

  (* norm is negative *)
  norm = -N[Flatten[{{0, 1}, list0}], prec];
  (* obj = N[Flatten[{{1, 0}, list0}], prec]; *)
  obj = ConstantArray[0, functionalCount];

  If[Length[norm] =!= functionalCount || Length[obj] =!= functionalCount,
    Print["Objective or normalization dimension mismatch."];
    Abort[]
  ];

  Print["functional dimension = ", functionalCount];
  Print["PMP blocks = ", Length[pols], " (expected ",
    2 + nPoints (Length[spins] + 1), ")"];
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
