test18Directory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{test18Directory, "..", "SDPB.m"}]];

m1 = N[2/5, 1000];
mA = N[1/5, 1000];
J1 = 0;
J2 = 2;
mgap = N[166/100, 1000];

(* Independent functional coordinates for the three null-sum-rule families. *)
crossNullIndices = Range[0, 4];
aaNullIndices = Range[0, 4];
bbNullIndices = Range[0, 4];

nullCount = Total[
  Length /@ {crossNullIndices, aaNullIndices, bbNullIndices}
];
list0 = ConstantArray[0, nullCount];


NABAB[n_, x_, J_, m_] := {(-m^2+x)^5/x^4,(J (7+J) (m^2-x)^4)/(4 x^3),(m^2-x)^4/x^4,-((J (7+J) (m^2-x)^3)/(4 x^3)),((-1+J) J (7+J) (8+J) (m^2-x)^2)/(40 x^2)}[[n+1]];


NAABB[n_, x_, J_, m_]:= {(J (7+J) x^(1/4) (-4 m^2+x)^(7/4) Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 Sqrt[x (-4 m^2+x)]),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(3/4),((-1+J) J (7+J) (8+J) (-4 m^2+x)^(3/4) Hypergeometric2F1[2-J,9+J,6,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(40 x^(3/4)),(J (7+J) (-4 m^2+x)^(3/4) Sqrt[x (-4 m^2+x)] Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 x^(7/4)),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(7/4)}[[n+1]];


NBBBB[n_, z_, J_, m_] := {z-1/2 J (7+J) z,1-1/20 J (7+J) (-13+J (7+J)),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z)),1/z^2-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^2),-((J (7+J) (-23+J (7+J)))/(20 z^2))}[[n+1]];


NAAAA[n_, x_, J_, mA_] := {((-4 mA^2+x)^(3/2) (32 mA^4+2 (-1+J) (8+J) mA^2 x-(-2+J (7+J)) x^2))/(2 x^(5/2)),1/(20 x^(7/2))Sqrt[-4 mA^2+x] (-1280 mA^6+960 mA^4 x+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 x^2-(-20+J (7+J) (-13+J (7+J))) x^3),1/(360 x^(9/2) Sqrt[-4 mA^2+x])(92160 mA^8-92160 mA^6 x+34560 mA^4 x^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 x^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) x^4),-((J (7+J) mA^2 (80 mA^4+2 (-38+J (7+J)) mA^2 x-(-23+J (7+J)) x^2))/(5 x^(7/2) Sqrt[-4 mA^2+x])),1/(10080 x^(11/2) (-4 mA^2+x)^(3/2))(-10321920 mA^10+12902400 mA^8 x-6451200 mA^6 x^2+1612800 mA^4 x^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 x^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) x^5)}[[n+1]];


NCross[n_, z_, J_] := {
  {0, 0, -1/2 NAABB[n, z, J, mA]},
  {0, NABAB[n, z, J, mA], 0},
  {-1/2 NAABB[n, z, J, mA], 0, 0}
};

NAAMatrix[n_, z_, J_] := {
  {NAAAA[n, z, J, mA], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

NBBMatrix[n_, z_, J_] := {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, NBBBB[n, z, J, mA]}
};

(* D=10 zero-subtracted positive moment for AA -> AA. *)
g0AAWeight[x_] :=
  (2 (-4 mA^2+x)^(7/2))/(Sqrt[x] (-2 mA^2+x));

(* In the block-diagonal basis {AA, AB, BB}, a standard universal spin-2
   couples equally to the two neutral channels and not to AB. *)
universalSpin2Direction = {1, 0, 1};

contractUniversalSpin2[matrix_] :=
  universalSpin2Direction . matrix . universalSpin2Direction;

neutralBlock[matrix_] := matrix[[{1, 3}, {1, 3}]];


(* ---------------------------------------------------------------------- *)
(* Large-spin diagnostics (not yet imposed as PMP constraints)             *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  phaseABAB10,
  aabbFixedMassGrowthBase10,
  fixedMassLargeJData10,
  largeJImpactWeights10
];

phaseABAB10[x_, m_] := (x - m^2)^7/x^4;


(* At fixed x > 4 mass^2, the AABB Gegenbauer polynomial is evaluated
   outside [-1,1] and grows exponentially with spin. *)
aabbFixedMassGrowthBase10[x_, mass_] :=
  (Sqrt[x] + 2 mass)/Sqrt[x - 4 mass^2];

fixedMassLargeJData10[x_, mass_] := <|
  "ABABCoefficientAfterJ4Scaling" -> 1/(40 x^2),
  "AABBGrowthBasePerSpin" -> aabbFixedMassGrowthBase10[x, mass],
  "AABBEvenSpinRatio" -> aabbFixedMassGrowthBase10[x, mass]^2
|>;

(* Joint large-J/large-mass limit with b = J/Sqrt[x] fixed.  If the
   physical convention is bPhysical = 2 J/Sqrt[x], use b -> bPhysical/2. *)
largeJImpactWeights10[b_, mass_] := <|
  "ABAB" -> b^4/40,
  "AABB" -> b^4 Hypergeometric0F1[6, mass^2 b^2]/40
|>;


LaunchKernels[];


PMP2SDP[datfile_, prec_:600] := Module[
    {
        npts, phiSamples, massSamples, Jmax,
        evenSpinSamples, oddSpinSamples,
        Poly, Poly1st, PolyABOdd, Poly2nd, pols, norm, obj,
        functionalCount, expectedBlocks
    },
    If[! TrueQ[Min[m1, 1, mgap] > 4 mA^2],
      Print["Invalid spectrum: every sampled pole must satisfy z > 4 mA^2."];
      Abort[]
    ];

    (* Paper eq. (D.1): Chebyshev nodes in the conformal angle.
       Here m = 1-mgap/z = Sin[phi/2]^2. *)
    npts = 200;
    phiSamples = N[Table[
      Pi/2 + Pi/2 Cos[(k + 1/2) Pi/npts],
      {k, 0, npts - 1}
    ], prec];
    massSamples = Sin[#/2]^2 & /@ phiSamples;

    If[
      Length[massSamples] =!= npts ||
        ! AllTrue[massSamples, 0 < # < 1 &],
      Print["Invalid mass-sampling grid."];
      Abort[]
    ];

    (* Neutral AA/BB states have even spin.  The independent AB spectral
       sector also admits odd spin.  Large-J constraints are not imposed yet. *)
    Jmax = 100;
    evenSpinSamples = Range[0, Jmax, 2];
    oddSpinSamples = Range[1, Jmax, 2];

    (* continuous spectrum *)
    Poly[j_, x_, y_] := Module[{g0, lambda22, polys},
      (* Block-diagonal basis {AA, AB, BB}: the neutral {AA, BB} sector is
         a 2 x 2 block, while AB is an independent 1 x 1 block. *)
      g0 = {{g0AAWeight[x], 0, 0}, {0, 0, 0}, {0, 0, 0}};

      (* Only the isolated second state contributes to lambda22. *)
      lambda22 = ConstantArray[0, {3, 3}];

      polys = Join[
        {g0, lambda22},
        NCross[#, x, j] & /@ crossNullIndices,
        NAAMatrix[#, x, j] & /@ aaNullIndices,
        NBBMatrix[#, x, j] & /@ bbNullIndices
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 3}, {column, 3}
        ]
       ]
    ];

    (* The first isolated state is neutral, so it only sees {AA, BB}. *)
    Poly1st[j_, x_, y_] := Module[{g0, lambda22, polys},
      g0 = {{g0AAWeight[x], 0, 0}, {0, 0, 0}, {0, 0, 0}};
      lambda22 = ConstantArray[0, {3, 3}];

      polys = neutralBlock /@ Join[
        {g0, lambda22},
        NCross[#, x, j] & /@ crossNullIndices,
        NAAMatrix[#, x, j] & /@ aaNullIndices,
        NBBMatrix[#, x, j] & /@ bbNullIndices
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
          Table[polys[[k, row, column]], {k, Length[polys]}],
          {row, 2}, {column, 2}
        ]
      ]
    ];

    (* Odd spins occur only in the independent AB spectral sector. *)
    PolyABOdd[j_, x_, y_] := Module[{g0, lambda22, polys},
      g0 = 0;
      lambda22 = 0;

      polys = Join[
        {g0, lambda22},
        NABAB[#, x, j, mA] & /@ crossNullIndices,
        ConstantArray[0, Length[aaNullIndices] + Length[bbNullIndices]]
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        {{polys}}
      ]
    ];

    (* The isolated spin-2 state has the fixed universal neutral coupling
       direction {gHAA, gHAB, gHBB} proportional to {1, 0, 1}.  Its only
       nonnegative spectral variable is the bare coupling squared lambdaH2. *)
    Poly2nd[j_, x_, y_] := Module[{g0, lambdaH2, polys},
      g0 = {{g0AAWeight[x], 0, 0}, {0, 0, 0}, {0, 0, 0}};

      (* Unit coefficient means that the second functional coordinate
         normalizes the bare universal coupling squared. *)
      lambdaH2 = 1;

      polys = Join[
        {contractUniversalSpin2[g0], lambdaH2},
        contractUniversalSpin2[NCross[#, x, j]] & /@ crossNullIndices,
        contractUniversalSpin2[NAAMatrix[#, x, j]] & /@ aaNullIndices,
        contractUniversalSpin2[NBBMatrix[#, x, j]] & /@ bbNullIndices
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        {{polys}}
      ]
    ];
    
    pols = Flatten[{
      Flatten[ N[ ParallelTable[ Poly1st[i, m1, x], {i, J1, J1, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Poly2nd[i, 1, x], {i, J2, J2, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Poly[i, mgap/(1-m), x], {i, evenSpinSamples}, {m, massSamples}], prec] ],
      Flatten[ N[ ParallelTable[ PolyABOdd[i, mgap/(1-m), x], {i, oddSpinSamples}, {m, massSamples}], prec] ]
    }, 1];

    expectedBlocks = 2 + npts (
      Length[evenSpinSamples] + Length[oddSpinSamples]
    );
    If[Length[pols] =!= expectedBlocks,
      Print[
        "Unexpected block count: built ", Length[pols],
        ", expected ", expectedBlocks, "."
      ];
      Abort[]
    ];

    Print[
      "Built ", Length[pols],
      " numerical PMP blocks: neutral and AB/even J = 0, 2, ..., ", Jmax,
      "; AB/odd J = 1, 3, ..., ", Jmax - 1, "."
    ];

    (* Bound the universal isolated-state coupling squared in units of
       the AA positive moment g0AA. *)
    norm = -1 * N[Flatten[{{0, 1}, list0}], prec];
    obj = -1 * N[Flatten[{{1, 0}, list0}], prec];

    functionalCount = 2 + Length[list0];
    If[Length[norm] =!= functionalCount || Length[obj] =!= functionalCount,
      Print["Objective/normalization dimension mismatch."];
      Abort[]
    ];

    Print["size of norm = ", Length[norm]];
    Print["size of obj = ", Length[obj]];
    Print["Writing ", datfile, "..."];
    WritePmpJson[datfile, SDP[obj, norm, pols], prec, getAnalyticSampleData];

    Print["Wrote ", datfile, "."]
];

PMP2SDP[FileNameJoin[{test18Directory, "n_pmp.json"}], 1000];
