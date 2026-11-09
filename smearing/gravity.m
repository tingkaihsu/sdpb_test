FileDir = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{FileDir, "..", "SDPB.m"}]];

(* partial-wave basis *)
d = 6;
calP[J_, z_] := Hypergeometric2F1[-J, J + d - 3, (d - 2)/2, (1 - z)/2];

(* improved sum rules *)
kappa[J_] := 
  Together[(2 + v)/(1 + v)^2 calP[J, 1 + 2 v] - 
    v^2 ((4 + 3 v)/(1 + v)^2 calP[J, 1] + (4 v)/(
        1 - v^2) (D[calP[J, z], z] /. {z -> 1}))];
(* null constraints *)
chi[k_, J_] := 
  Together[
   Module[{w}, (2 + 
        v) calP[J, 1 + 2 v]/(v (1 + v))/(v (1 + v))^(k/2) - 
     SeriesCoefficient[(2 + w) (1 - 
         w) (1 + 2 w)/((v - w) w (1 - v) (1 + w) (1 + v + w))*
       calP[J, 1 + 2 w]/(w (1 + w))^(k/2), {w, 0, -1}]]];


QR[expr_] := Module[{num = Numerator[expr], den = Denominator[expr], R}, 
    If[FreeQ[den, v], {CoefficientList[expr, v], 0}, 
        R = Together[(1 - v) expr] /. v -> 1;
        {CoefficientList[PolynomialQuotient[num - R den/(1 - v), den, v], v], R}
    ]
];

smearQR[{qc_, residue_}, n_, q_] := (
  Sum[qc[[j + 1]] (-1)^j q^(-j)/(n + 2 j + 1),
      {j, 0, Length[qc] - 1}] +
  residue Hypergeometric2F1[1, (n + 1)/2, (n + 3)/2, -1/q]/(n + 1)
);

(* smeared, improved sum rule *)
rowF[J_, m2_, fList_] := Table[smearQR[QR[kappa[J]], n, m2], {n, fList}]/m2^2;
(* smeared null constraints *)
rowX[k_, J_, m2_, iList_] := Table[smearQR[QR[chi[k, J]], i, m2], {i, iList}] m2^(-1 - 3 k/2);


LaunchKernels[];


g2Fixed = 0;

PMP2SDP[datfile_, prec_: 1000] := Module[
    {
        npts = 150, stateBlock, bBlock, fList, bigB, mmax, nNull, masslist,
        largeBBlock, Acoef, trunc, Nexp, ent, osc, nullZeros, cG, c2, c3,
        phiSamples, massSamples, massList, norm, obj
    },


    chebMasses[np_, prec_] := Module[{phi = Table[Pi/2 + Pi/2 Cos[(k + 1/2) Pi/np], {k, 0, np -1}]},
        N[1/(1 - Sin[#/2]^2) & /@ phi, prec]
    ];
    suppMasses[J_, dphi_, bmax_, prec_] := Module[{nn = J + 3/2, K},
        K = Floor[N[Min[nn Pi, bmax]/dphi, 50]];
        Table[N[1/Sin[j dphi/(2 nn)]^2, prec], {j, 1, K}]
    ];

    

    fList = Range[3/2, 15/2]; bigB = 40; mmax = 2; nNull = 10;

    massList[J_] := Join[
        chebMasses[npts, prec],
        suppMasses[J, dphi, bmax, prec]
    ];
    Acoef[n_] := 2^n Gamma[(n + 1)/2]/Gamma[(3 - n)/2];
    trunc[e_] := Normal[Series[e bb^(5/2), {bb, Infinity, mmax}]]/bb^(5/2);
    Nexp = Max[fList] + 1;

    osc[n_, mmax_] := Module[{ser, c1, c2, al, be},
        ser = Expand[ Normal[ Series[ HypergeometricPFQ[{(n + 1)/2}, {2, (n + 3)/2}, -bb^2/4]/(n + 1),
            {bb, Infinity, mmax + 3} ] ] ] /. E^(y_) :> E^Expand[y];

        c1 = Coefficient[ser, E^(I bb)]; c2 = Coefficient[ser, E^(-I bb)];
        al = c1 E^(5 I Pi/4); be = c2 E^(-5 I Pi/4);
        {FullSimplify[ExpandAll[al + be]], FullSimplify[ExpandAll[I (al - be)]]}
    ];

    ent[n_] := Module[{Bn, Cn}, 
        {Bn, Cn} = trunc /@ osc[n, mmax];
        Expand[bb^Nexp {Acoef[n] bb^(-n - 1) + Bn, Cn, Acoef[n] bb^(-n - 1) - Bn} /. bb -> bigB + x]
    ];


    stateBlock[J_, m2_] := PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, x],
        {{Join[ rowF[J, m2, fList], rowX[4, J, m2, Range[0, 5] ], rowX[6, J, m2, Range[0, 3] ] ]}}
    ];
    
    bBlock[b_] := PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, x],
        {{Join[ Table[ HypergeometricPFQ[ {(n + 1)/2}, {2, (n + 3)/2}, -b^2/4]/(n + 1), {n, fList} ], ConstantArray[0, nNull] ]}}
    ];
    largeBBlock := Module[{e = ent /@ fList, z = ConstantArray[0, nNull]},
        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, x],
            {{Join[e[[All, 1]], z], Join[e[[All, 2]], z]},
            {Join[e[[All, 2]], z], Join[e[[All, 3]], z]}}
        ]
    ];


    pols = Join[
        Flatten[ParallelTable[stateBlock[J, m2], {J, 0, 42, 2}, {m2, massList[J]}]],
        Table[bBlock[b], {b, Range[1/250, 40, 1/32]}],
        {largeBBlock}
    ];
    
    nullZeros = ConstantArray[0, nNull];

    cG = Join[Table[1/(n - 1), {n, fList}], nullZeros];
    c2 = Join[Table[2/(n + 1), {n, fList}], nullZeros];
    c3 = Join[Table[1/(n + 3), {n, fList}], nullZeros];
    norm = -c3;
    obj = -(cG + g2Fixed c2);
    WritePmpJson[datfile, SDP[obj, norm, N[pols, prec]], prec, getAnalyticSampleData];
];    


PMP2SDP["n_pmp.json", 1000]