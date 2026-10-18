(* ::Package:: *)

(* ======================================================================== *)
(* mixed_LT_lambdaAB.m                                                      *)
(*                                                                          *)
(* Linear-trajectory bootstrap of the mixed A (massive, Z2 +) / B           *)
(* (massless, Z2 -) system, D = 10.  Bound on                               *)
(*                                                                          *)
(*      lambda_AA lambda_BB / Sqrt[ghatA ghatB]   at fixed r = ghatA/ghatB  *)
(*                                                                          *)
(* for the isolated X2 (m2sq, J2).  Notation: mixed_state_bootstrap_v2.tex. *)
(*                                                                          *)
(*   v   = {ghatB, ghatA, lambda22, nAAAA, nBBBB, nAB}                      *)
(*   a   = (a1, a2, -1, alpha, beta, gamma)          (upper bound)          *)
(*   a.v >= 0 (PSD) on every state,  maximize a.obj = -(a1 + a2 r)          *)
(*   SDPB objective V  =>  lambda lambda/Sqrt[ghatA ghatB] <= -V/Sqrt[r]    *)
(*                                                                          *)
(* Usage:  wolframscript -file mixed_LT_lambdaAB.m 1/5 1                    *)
(*         (mAsq, r; exact numbers.  Defaults 1/5 and 1.)                   *)
(* ======================================================================== *)

LTDirectory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{LTDirectory, "..", "SDPB.m"}]];


(* ------------------------------------------------------------------------ *)
(* BLOCK 1.  Input                                                          *)
(* ------------------------------------------------------------------------ *)

mAsq = If[Length[$ScriptCommandLine] >= 2,
  ToExpression[StringTrim[$ScriptCommandLine[[2]]]], 1/5];
rRatio = If[Length[$ScriptCommandLine] >= 3,
  ToExpression[StringTrim[$ScriptCommandLine[[3]]]], 1];

mA     = Sqrt[mAsq];
JA     = 0;
m2sq   = 1;        (* X2, units m2 = 1                                    *)
J2     = 2;
mgapSq = 2;        (* neutral continuum M^2 >= mgapSq, J even             *)
mcSq   = 1;        (* charged continuum M^2 >= mcSq,   all J              *)

boundType = "upper";   (* "upper" or "lower"                              *)

If[!(MatchQ[mAsq, _Integer | _Rational] && MatchQ[rRatio, _Integer | _Rational] &&
     0 < mAsq && 4 mAsq < m2sq < mgapSq && mAsq < mcSq && rRatio > 0),
  Print["Invalid input: need exact 0 < 4 mAsq < m2sq < mgapSq, mAsq < mcSq, r > 0."];
  Exit[2]
];
Print["mAsq = ", mAsq, ", r = ", rRatio, ", bound = ", boundType];


(* ------------------------------------------------------------------------ *)
(* BLOCK 2.  Null indices                                                   *)
(*                                                                          *)
(* row n  <->  (k,q), the coefficient of (s-s0)^(k-q) t^q                   *)
(*   AAAA  (s0 = 0)      n = 0..3 : (1,0) (2,0) (3,0) (3,1)  = test20       *)
(*   BBBB  (s0 = 0)      n = 0..2 : (1,0) (2,0) (3,0)                       *)
(*         [(3,1) vanishes identically for massless externals]              *)
(*   AB    (s0 = -mA^2)  n = 0..5 : (0,0) (1,0) (1,1) (2,0) (2,1) (2,2)     *)
(* ------------------------------------------------------------------------ *)
g0AAWeight[x_] := 2/(x - 2 mAsq);     (* ghatA at (s,t) = (2mA^2, 0) *)
g0BBWeight[x_] := 2/x;                (* ghatB at (s,t) = (0, 0)     *)

(* phi = (2p)^7/Sqrt[x].  Neutral element ij carries Sqrt[phi_i phi_j]:   *)
(* congruence diag(Sqrt[phiAA], Sqrt[phiBB]), PSD cone unchanged.         *)
phiAA[x_] := (x - 4 mAsq)^(7/2)/Sqrt[x];
phiBB[x_] := x^3;
phiAB[x_] := (x - mAsq)^7/x^4;

aaNullIndices = Range[0, 3];
bbNullIndices = Range[0, 2];
abNullIndices = Range[0, 1];

nullCount = Total[
  Length /@ {abNullIndices, aaNullIndices, bbNullIndices}
];

list0 = ConstantArray[0, nullCount];
functionalDimension = 3 + nullCount;


NABAB[n_,x_,J_,m_] := {}[[n+1]];


NAABB[n_,x_,J_,m_] := {}[[n+1]];


NBBBB[n_,z_,J_,m_] :={}[[n+1]];


NAAAA[n_,x_,J_,mA_] := {}[[n+1]];



NABMatrix[n_, z_, J_] := {
  {0, -1/2 NAABB[n, z, J, mA]},
  {-1/2 NAABB[n, z, J, mA], 0}
};

NAAMatrix[n_, z_, J_] := {
  {NAAAA[n, z, J, mA], 0},
  {0,0}
};

NBBMatrix[n_, z_, J_] := {
  {0, 0},
  {0, NBBBB[n, z, J, mA]}
};



(* light pole terms *)

nlABAB[n_,x_,J_,m_] := {}[[n+1]];



nlAABB[n_,x_,J_,m_] := {}[[n+1]];



nlBBBB[n_,z_,J_,m_] := {}[[n+1]];



nlAAAA[n_,x_,J_,mA_] := {}[[n+1]];



(* light pole terms *)

nlABMatrix[n_, z_, J_] := {
  {0, -1/2 nlAABB[n, z, J, mA]},
  {-1/2 nlAABB[n, z, J, mA], 0}
};


nlAAMatrix[n_, z_, J_] := {
  {nlAAAA[n, z, J, mA], 0},
  {0,0}
};


nlBBMatrix[n_, z_, J_] := {
  {0, 0},
  {0, nlBBBB[n, z, J, mA]}
};




(* ------------------------------------------------------------------------ *)
(* BLOCK 12.  Assemble and write                                            *)
(* ------------------------------------------------------------------------ *)

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
  {
    npts = 150, Jmax = 100, sgn,
    Fn, Fc, Fn2, FA, pols,
    phiSamples, massSamples, evenSpins, allSpins, norm, obj
  },

  phiSamples = N[Table[Pi/2 + Pi/2 Cos[(k + 1/2) Pi/npts], {k, 0, npts - 1}], prec];
  massSamples = Sin[#/2]^2 & /@ phiSamples;
  evenSpins = Range[0, Jmax, 2];
  allSpins  = Range[0, Jmax];

  Fn[j_, x_, y_] := Module[{w, gB, gA, lam, polys},
    w  = Sqrt[phiAA[x]*phiBB[x]];                     (* off-diagonal weight *)
    gB = {{0, 0}, {0, phiBB[x]*g0BBWeight[x]}};
    gA = {{phiAA[x]*g0AAWeight[x], 0}, {0, 0}};

    lam = {{0, 0}, {0, 0}};   (* lam0 = {{0,1/2},{1/2,0}} for X2, zero matrix otherwise *)

    polys = Join[
      {gB, gA, lam},
      w*NABMatrix[#, x, j] & /@ abNullIndices,
      phiAA[x]*NAAMatrix[#, x, j] & /@ aaNullIndices,
      phiBB[x]*NBBMatrix[#, x, j] & /@ bbNullIndices
    ];

    PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 2}, {column, 2}
        ]
    ]
  ];

  Fc[j_, x_, y_] := Module[{polys},
    polys = Join[
      {0, 0, 0},
      phiAB[x]*NABAB[#, x, j, mA] & /@ abNullIndices,
      0 & /@ aaNullIndices,
      0 & /@ bbNullIndices
    ];

    PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        {{polys}}
    ]
  ];

  (* second isolated state *)
  Fn2[j_, x_, y_] := Module[{w, gB, gA, lam, polys},
    w  = Sqrt[phiAA[x]*phiBB[x]];                     (* off-diagonal weight *)
    gB = {{0, 0}, {0, phiBB[x]*g0BBWeight[x]}};
    gA = {{phiAA[x]*g0AAWeight[x], 0}, {0, 0}};

    lam = w*{{0, 1/2}, {1/2, 0}};   (* lam0 = {{0,1/2},{1/2,0}} for X2, zero matrix otherwise *)

    polys = Join[
      {gB, gA, lam},
      w*NABMatrix[#, x, j] & /@ abNullIndices,
      phiAA[x]*NAAMatrix[#, x, j] & /@ aaNullIndices,
      phiBB[x]*NBBMatrix[#, x, j] & /@ bbNullIndices
    ];

    PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 2}, {column, 2}
        ]
    ]
  ];

  (* light A and B, acting on (gAAA, gABB); B shares gABB -> (2,2) slot *)
  FA[j_, x_, y_] := Module[{zero2, polys},
    zero2 = ConstantArray[0, {2, 2}];

    polys = Join[
      {zero2, zero2, zero2},
      nlABMatrix[#, mAsq, JA] + {{0, 0}, {0, nlABAB[#, 0, 0, mA]}} & /@ abNullIndices,
      nlAAMatrix[#, mAsq, JA] & /@ aaNullIndices,
      nlBBMatrix[#, mAsq, JA] & /@ bbNullIndices
    ];

    PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 2}, {column, 2}
        ]
    ]
  ];


  pols = Flatten[{
      Flatten[ N[ ParallelTable[ FA[i, mA, x], {i, JA, JA, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Fn2[i, 1, x], {i, J2, J2, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Fn[i, mgapSq/(1-m), x], {i, evenSpins}, {m, massSamples}], prec] ],
      Flatten[ N[ ParallelTable[ Fc[i, 1/(1-m), x], {i, allSpins}, {m, massSamples}], prec] ]
  }, 1];

  
  sgn  = If[boundType === "upper", 1, -1];
  norm = -sgn N[Join[{0, 0, 1}, list0], prec];
  obj  = -N[Join[{1, rRatio, 0}, list0], prec];

  If[
    Length[norm] =!= functionalDimension ||
    Length[obj]  =!= functionalDimension,
    Print["Objective or normalization dimension mismatch."];
    Abort[]
  ];

  Print["functional dimension = ", functionalDimension,
        "  (3 + ", Length[abNullIndices], " AB + ", Length[aaNullIndices],
        " AAAA + ", Length[bbNullIndices], " BBBB)"];
  Print["number of PMP blocks = ", Length[pols]];
  Print["bound = ", If[sgn == 1, "-V/Sqrt[r]", "+V/Sqrt[r]"],
        " with V the SDPB objective; Sqrt[r] = ", N[Sqrt[rRatio], 10]];

  WritePmpJson[
    datfile,
    SDP[obj, norm, pols],
    prec,
    getAnalyticSampleData
  ]
];

PMP2SDP[FileNameJoin[{LTDirectory,"n_pmp.json"}], 1000];
