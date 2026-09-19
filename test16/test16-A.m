(* ::Package:: *)

sourceDirectory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{sourceDirectory, "..", "SDPB.m"}]];

(* Keep the spectral data exact until the final numerical conversion. *)
m1 = 1/2;
mgap = 83/50;
mA = 1/5;

J1 = 0;
J2 = 2;

NullIndices = Range[0, 4];
nulllist = {4, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];
functionalDimension = 2 + Length[list0];


Nlist[n_, x_, J_, mA_] := {1/2 x^3 (-4 mA^2+x)^3 (-2 mA^2+x) (32 mA^4+2 (-1+J) (8+J) mA^2 x-(-2+J (7+J)) x^2),1/20 x^2 (-4 mA^2+x)^2 (-2 mA^2+x) (-1280 mA^6+960 mA^4 x+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 x^2-(-20+J (7+J) (-13+J (7+J))) x^3),1/360 x (-4 mA^2+x) (-2 mA^2+x) (92160 mA^8-92160 mA^6 x+34560 mA^4 x^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 x^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) x^4),-(1/5) J (7+J) mA^2 x^2 (-4 mA^2+x) (-2 mA^2+x) (80 mA^4+2 (-38+J (7+J)) mA^2 x-(-23+J (7+J)) x^2),1/10080(-2 mA^2+x) (-10321920 mA^10+12902400 mA^8 x-6451200 mA^6 x^2+1612800 mA^4 x^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 x^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) x^5)}[[n+1]];


(* ---------------------------------------------------------------------- *)
(* PMP assembly                                                           *)
(* ---------------------------------------------------------------------- *)

(* Dense low spins locate the usual extremal states.  The sparse probes
   monitor the transition to the analytic asymptotic blocks.  A final
   solution must still be checked for omitted finite-spin violations. *)
Jmax = 200;
JList = {250, 300, 400, 500, 700, 1000, 1500, 2500, 5000};
continuumJList = DeleteDuplicates @ Join[
  Range[0, Jmax, 2],
  JList
];

(* Positive blockwise rescaling for the q^4 = O(J^8) growth, where
   q = J (J + 7) is the ten-dimensional spin Casimir. *)
spinScale[J_Integer] := 1/(1 + J (J + 7))^4;

(* Analytic boundaries of the large-spin spectrum.  The fixed-mass
   vector was divided by its strictly positive mass-dependent factor. *)
fixedMassLargeJVector = {0, 0, 0, 0, 0, 0, -1};

impactVector[r_] := {
  2,
  0,
  -r/2,
  -r^2/20,
  -r^3/360,
  0,
  -r^4/10080
};

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
    {
      continuumBlocks, specialBlocks, asymptoticBlocks,
      pols, norm, obj, expectedBlocks
    },

    Poly[j_, x_, y_] := Module[{g0, lambda22, polys},
        g0 = 2 x^5 (-4 mA^2 + x)^5;

        (* Only the isolated second state contributes to lambda22. *)
        lambda22 = 0;

        polys = Join[
            {g0, lambda22},
            Nlist[#, x, j, mA] & /@ NullIndices
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{spinScale[j] polys}}
        ]
    ];

    (* The first isolated state. *)
    Poly1st[j_, x_, y_] := Module[{g0, lambda22, polys},
        g0 = 2 x^5 (-4 mA^2 + x)^5;
        lambda22 = 0;

        polys = Join[
            {g0, lambda22},
            Nlist[#, x, j, mA] & /@ NullIndices
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{polys}}
        ]
    ];


    (* The isolated spin-2 state. *)
    Poly2nd[j_, x_, y_] := Module[{g0, lambdaH2, polys},
        g0 = 2 x^5 (-4 mA^2 + x)^5;
        lambdaH2 = 1;

        polys = Join[
            {g0, lambdaH2},
            Nlist[#, x, j, mA] & /@ NullIndices
        ];

        PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
            {{polys}}
        ]
    ];

    (* Fixed mass with J -> Infinity. *)
    FixedMassLargeJBlock[y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      {{fixedMassLargeJVector}}
    ];

    (* Joint limit with r = J (J + 7)/s fixed. *)
    ImpactBlock[y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      {{impactVector[y]}}
    ];

    specialBlocks = {
      Poly1st[J1, m1, x],
      Poly2nd[J2, 1, x]
    };

    (* One polynomial block per spin covers the complete mass interval
       s = mgap + x, x >= 0; no mass sampling is required. *)
    continuumBlocks = ParallelMap[
      Poly[#, mgap + x, x] &,
      continuumJList
    ];

    asymptoticBlocks = {
      FixedMassLargeJBlock[x],
      ImpactBlock[x]
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
    Print["number of PMP blocks = ", Length[pols] ];
    Print["finite J values = ", continuumJList];
    Print["Included fixed-mass and correlated large-J blocks."];

    WritePmpJson[
        datfile,
        SDP[obj, norm, pols],
        prec,
        getAnalyticSampleData
    ]
];

outputFile = FileNameJoin[{sourceDirectory, "n_pmp.json"}];
If[!TrueQ[$Test16SkipExport],
  PMP2SDP[outputFile, 600]
];
