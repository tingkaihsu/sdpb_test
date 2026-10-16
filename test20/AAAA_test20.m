(* ::Package:: *)

test20Directory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{test20Directory, "..", "SDPB.m"}]];

mAsq = 1/5;
mA = Sqrt[mAsq];
JA = 0;
m2sq = 1;
J2 = 2;
mgap = 2;


aaNullIndices = Range[0, 3];

nAA = Length[aaNullIndices];

list0 = ConstantArray[0, nAA];
functionalDimension = 2 + nAA;   (* {ghat, lambda22, nulls} *)
zA = ConstantArray[0, nAA];

(* poles q = 0 *)
KX[n_, x_, J_] := {(32*mA^4 - 16*mA^2*x + 14*J*mA^2*x + 2*J^2*mA^2*x + 2*x^2 - 7*J*x^2 - J^2*x^2)/(2*x^2*(-4*mA^2 + x)^2), (-1280*mA^6 + 960*mA^4*x - 240*mA^2*x^2 - 112*J*mA^2*x^2 + 82*J^2*mA^2*x^2 + 28*J^3*mA^2*x^2 + 2*J^4*mA^2*x^2 + 20*x^3 + 91*J*x^3 - 36*J^2*x^3 - 14*J^3*x^3 - J^4*x^3)/(20*x^3*(-4*mA^2 + x)^3), (92160*mA^8 - 92160*mA^6*x + 34560*mA^4*x^2 - 5760*mA^2*x^3 + 2016*J*mA^2*x^3 - 2260*J^2*mA^2*x^3 - 42*J^3*mA^2*x^3 + 242*J^4*mA^2*x^3 + 42*J^5*mA^2*x^3 + 2*J^6*mA^2*x^3 + 360*x^4 - 2142*J*x^4 + 1409*J^2*x^4 + 147*J^3*x^4 - 112*J^4*x^4 - 21*J^5*x^4 - J^6*x^4)/(360*x^4*(-4*mA^2 + x)^4), -1/5*(J*(7 + J)*mA^2*(80*mA^4 - 76*mA^2*x + 14*J*mA^2*x + 2*J^2*mA^2*x + 23*x^2 - 7*J*x^2 - J^2*x^2))/((4*mA^2 - x)^4*x^3)}[[n+1]];

NA[n_, x_, J_] := {1/mA^4, 1/mA^6, 1/mA^8, 0}[[n+1]];
(* D=10, AA -> AA at s = 2mA^2, t = 0, coupling normalization *)
g0AAWeight[x_] := 2/(-2*mA^2 + x);


(* ------------------------------------------------------------------------ *)
(* PHASE-SPACE WEIGHTS (heavy states only), D = 10: phi = (2p)^7/Sqrt[x]    *)
(*                                                                          *)
(* Element ij of a heavy block carries phi_ij on ALL its components         *)
(* (g0AA, g0BB, XH and every null row).  The off-diagonal weight is the     *)
(* geometric mean, so the PSD cone is unchanged (positive diagonal          *)
(* congruence D = diag(Sqrt[phiAA], Sqrt[phiBB])).                          *)
(*   phiBB g0BBWeight = 2 x^2  (the original BB forward weight)             *)
(*   phiX            = phiAABB of ABAB_check.m                              *)
(*   phiAB           = wAB of ABAB_check.m                                  *)
(* ------------------------------------------------------------------------ *)

phiAA[x_] := (x - 4 mA^2)^(7/2)/Sqrt[x];

spinScale[j_] := 1;   (* any positive constant per block; 1 = no rescaling *)


LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
    {
      npts, phiSamples, massSamples, pols, norm, obj
    },

    npts = 150;
    phiSamples = N[Table[
      Pi/2 + Pi/2 Cos[(k + 1/2) Pi/npts],
      {k, 0, npts - 1}
    ], prec];
    massSamples = Sin[#/2]^2 & /@ phiSamples;
    
    Jmax = 100;
    evenSpinSamples = Range[0, Jmax, 2];


    PolyX[j_, x_, y_] := Module[{g0, lam, polys},
        g0 = g0AAWeight[x];

        (* Only the isolated second state contributes to lambda22. *)
        lam = 0;

        polys = phiAA[x]*Join[
            {g0, lam},
            Table[KX[n,x,j], {n, aaNullIndices}]
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{spinScale[j] polys}}
        ]
    ];

    PolyA[j_, x_, y_] := Module[{g0, lam, polys},
        g0 = 0;   (* ghat = g00 + gAAA^2/mA^2 receives nothing from A *)

        (* Only the isolated second state contributes to lambda22. *)
        lam = 0;

        polys = Join[
            {g0, lam},
            Table[NA[n,x,j], {n, aaNullIndices}]
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{spinScale[j] polys}}
        ]
    ];

    PolyH[j_, x_, y_] := Module[{g0, lam, polys},
        g0 = g0AAWeight[x];

        (* Only the isolated second state contributes to lambda22. *)
        lam = 1;

        polys = phiAA[x]*Join[
            {g0, lam},
            Table[KX[n,x,j], {n, aaNullIndices}]
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{spinScale[j] polys}}
        ]
    ];

    pols = Flatten[{
      Flatten[ N[ ParallelTable[ PolyA[i, mAsq, x], {i, JA, JA, 2}], prec] ],
      Flatten[ N[ ParallelTable[ PolyH[i, m2sq, x], {i, J2, J2, 2}], prec] ],
      Flatten[ N[ ParallelTable[ PolyX[i, mgap/(1-m), x], {i, evenSpinSamples}, {m, massSamples}], prec] ]
    }, 1];


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
    Print["finite J values = ", evenSpinSamples];

    WritePmpJson[
        datfile,
        SDP[obj, norm, pols],
        prec,
        getAnalyticSampleData
    ]
];

PMP2SDP[FileNameJoin[{test20Directory, "n_pmp.json"}], 1000];