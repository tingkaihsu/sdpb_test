Import["../SDPB.m"];
m1 = N[2/5, 1000];
J1 = 0;
J2 = 2;
mgap = N[166/100, 1000];

(* 7 null constraints *)
nulllist = {6, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];

NBBBB[n_, z_, J_] := {z-1/2 J (7+J) z,1-1/20 J (7+J) (-13+J (7+J)),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z)),1/z^2-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^2),-((J (7+J) (-23+J (7+J)))/(20 z^2)),1/z^3-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)))/(403200 z^3),-((J (7+J) (540+J (7+J) (-53+J (7+J))))/(360 z^3))}[[n+1]];


NAAAA[n_, x_, J_]:= {((-4 mA^2+x)^(3/2) (32 mA^4+2 (-1+J) (8+J) mA^2 x-(-2+J (7+J)) x^2))/(2 x^(5/2)),1/(20 x^(7/2))Sqrt[-4 mA^2+x] (-1280 mA^6+960 mA^4 x+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 x^2-(-20+J (7+J) (-13+J (7+J))) x^3),1/(360 x^(9/2) Sqrt[-4 mA^2+x])(92160 mA^8-92160 mA^6 x+34560 mA^4 x^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 x^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) x^4),-((J (7+J) mA^2 (80 mA^4+2 (-38+J (7+J)) mA^2 x-(-23+J (7+J)) x^2))/(5 x^(7/2) Sqrt[-4 mA^2+x])),1/(10080 x^(11/2) (-4 mA^2+x)^(3/2))(-10321920 mA^10+12902400 mA^8 x-6451200 mA^6 x^2+1612800 mA^4 x^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 x^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) x^5),1/(180 x^(9/2) (-4 mA^2+x)^(3/2))J (7+J) (11520 mA^8-11520 mA^6 x-4 (-936+J (7+J) (-26+J (7+J))) mA^4 x^2+2 (-216+J (7+J) (-26+J (7+J))) mA^2 x^3-9 (-23+J (7+J)) x^4),1/(403200 x^(13/2) (-4 mA^2+x)^(5/2))(1651507200 mA^12-2477260800 mA^10 x+1548288000 mA^8 x^2-516096000 mA^6 x^3+96768000 mA^4 x^4+2 (-4838400+(-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J)) mA^2 x^5-(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) x^6)}[[n+1]];


Nlist[n_,z_,J_] := {
	{NAAAA[n,z,J],0,0},
	{0,NBBBB[n,z,J],0},
	{0,0,0}
};


polyify[expr_] := Expand @ Cancel @ Together[expr];

PolyInfBBBB[n_,J_,x_] := {0,0,0,0,0,-(1/(403200 z^3)),0}[[n+1]];


PolyInfAAAA[n_,J_,x_] := {0,0,0,0,0,0,(2 mA^2-x)/(403200 x^(3/2) (-4 mA^2+x)^(5/2))}[[n+1]];


NPolyInf[n_,J_,x_] := {
  {PolyInfAAAA[n,J,x],0,0},
  {0,PolyInfBBBB[n,J,x],0},
  {0,0,0}
};

LaunchKernels[];


PMP2SDP[datfile_, prec_:600] := Module[
    {
        xTiers, jTiers, Poly, PolyInf,
        Poly2nd, pols, norm, obj,
        functionalCount, functionalCovered, missingFunctionals
    },

    xTiers = {
    10^Subdivide[-4, -3, 49],
    10^Subdivide[-3, -2, 49],
    10^Subdivide[-2, -1, 49],
    10^Subdivide[-1, Log10[1 - 10^-4], 49]
    };

    nPerTier = 50;
    n = 4 nPerTier;
    jMax = 50000;
    alp = 3/2;

    jAll = DeleteDuplicates[
      Round[jMax (Range[0, n - 1]/(n - 1))^alp]
    ];

    jAll = Sort[jAll];

    jTiers = Partition[jAll, nPerTier];

    (* continuous spectrum *)
    Poly[j_, x_, y_] := Module[{g0, lambda22, polys},
      (* normalization constant on BBBB sector *)
      g0 = {{0, 0, 0}, {0, x^3*2/x, 0}, {0, 0, 0}};

      (* state 2 on-shell couplings *)
      lambda22 = {{(-4 mA^2+x)^(7/2)/Sqrt[x]*0, 0, 0}, {0, x^3*0, 0}, {0, 0, 0}};

      polys = Join[
        {g0, lambda22},
        Table[ Nlist[n, x, j] , {n, 0, nulllist[[1]]} ]
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 3}, {column, 3}
        ]
       ]
    ];

    PolyInf[j_, x_, y_] := Module[{polys, n0},
      n0 = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      polys = Join[
        {n0, n0},
        Table[ NPolyInf[n, j, x], {n, 0, nulllist[[1]]} ]
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 3}, {column, 3}
        ]
      ]
    ];

    (* poly of the second state *)
    Poly2nd[j_, z_, y_] := Module[{g0, lambda22, polys},
      
      g0 = {{0, 0, 0}, {0, z^3*2/z, 0}, {0, 0, 0}};
      lambda22 = {{(-4 mA^2+z)^(7/2)/Sqrt[z]*1, 0, 0}, {0, z^3*1, 0}, {0, 0, 0}};
      polys = Join[
        {g0, lambda22},
        Table[ Nlist[n, z, j],
          {n, 0, nulllist[[1]]}
        ] 
      ];
      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 3}, {column, 3}
        ]
      ]
    ];
    
    pols = Flatten[{
      Flatten[ N[ ParallelTable[ Poly[i, m1, x], {i, J1, J1, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Poly2nd[i, 1, x], {i, J2, J2, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Poly[i, mgap*1/(1-m), x], {i, Flatten[jTiers]}, {m, Flatten[xTiers]}], prec] ],
      Flatten[ N[ ParallelTable[ PolyInf[i, mgap*1/(1-m), x], {i, 0, 0, 2}, {m, Flatten[xTiers]}], prec] ]
    }, 1];

    Print["Built ", Length[pols], " numerical PMP blocks."];

    norm = -1 * N[Flatten[{{0, 1}, list0}], prec];
    obj = -1 * N[Flatten[{{1, 0}, list0}], prec];

    Print["size of norm = ", Length[norm]];
    Print["size of obj = ", Length[obj]];
    Print["Writing ", datfile, "..."];
    WritePmpJson[datfile, SDP[obj, norm, pols], prec, getAnalyticSampleData];

    Print["Wrote ", datfile, "."]
];

PMP2SDP["n_pmp.json", 1000];
