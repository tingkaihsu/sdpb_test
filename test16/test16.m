(* ::Package:: *)

Import["../SDPB.m"]
m1 = N[45/100, 1000];
m2 = N[1, 1000];
J1 = 0;
mgap = N[166/100, 1000];

nulllist = {8, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];

Nlist[n_, x_, J_] := {(2-J (7+J))/z^2,1/z^3+((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) m1^2 m2^2 (7 m1^2-16 m2^2)-36 J (7+J) (-13+J (7+J)) (m1^4+8 m2^4) z-180 (-2+J (7+J)) (2 m1^2+7 m2^2) z^2)/(720 (m1^4+8 m2^4) z^4),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(180 z^4)),(-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))) m1^2 m2^2 (m1^4+8 m2^4))-28 (-2880 m2^6+(-2+J) J (7+J) (9+J) (-17+J (7+J)) (m1^6+8 m2^6)) z+5040 (-2+J (7+J)) (m1-m2) (m1+m2) z^3+10080 m1^2 (8 m2^6+m1^4 (m2^2+z)))/(10080 m1^2 m2^2 (m1^4+8 m2^4) z^5),1/(80 z^5) (-4 J (7+J) (-23+J (7+J))+1/(m2^2 (m1^4+8 m2^4)) z (-((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) m1^4)+180 (-2+J (7+J)) z^2)),-(1/(403200 z^6))(-403200+3769920 J-3778144 J^2+504980 J^3+368780 J^4-41405 J^5-18207 J^6+70 J^7+370 J^8+35 J^9+J^10+(1120 (-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) (m1^8+8 m2^8) z^2)/(m1^4 m2^4 (m1^4+8 m2^4))-(201600 (-2+J (7+J)) (m1^4-m2^4) z^4)/(m1^4 m2^4 (m1^4+8 m2^4))),1/(720 z^6) (-2 J (7+J) (540+J (7+J) (-53+J (7+J)))+1/(m2^4 (m1^4+8 m2^4)) 9 z^2 (-((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) m1^4)+180 (-2+J (7+J)) z^2)),1/(720 z^6) (2 J (7+J) (540+J (7+J) (-53+J (7+J)))+1/(m2^4 (m1^4+8 m2^4)) (9 (-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) m1^4 z^2-1620 (-2+J (7+J)) z^4)),1/z^7-1/(21772800 z^7) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))-((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) (m1^10+8 m2^10))/(360 m1^6 m2^6 (m1^4+8 m2^4) z^4)+((-2+J (7+J)) (m1^6-m2^6))/(2 m1^6 m2^6 (m1^4+8 m2^4) z^2)}[[n+1]];

polyify[expr_, var_] := Expand @ Cancel @ Together[expr];

Poly[J_, z_, y_] := Module[{pref, polys, first},
  pref = z^7;

  first = polyify[pref*-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) m1^4 m2^4 (64 m1^2-19 m2^2)+180 (-2+J (7+J)) (8 m1^6+19 m2^6) z^2-5760 (m1^4+8 m2^4) z^3)/(2880 (m1^4+8 m2^4) z^4)), z];
  second = polyify[pref*(m2^6 (-((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) m1^4)+180 (-2+J (7+J)) z^2))/(360 (m1^4+8 m2^4) z^4), z]
  

  polys = Table[
    polyify[pref*Nlist[n, z, J], z]
    , {n, 0, nulllist[[1]]}
  ];

  If[!AllTrue[Join[{first, second}, polys], PolynomialQ[#, z] &],
    Print["Non-polynomial terms remain."];
    Print[Pick[Range[1 + Length[polys]], Not /@ (PolynomialQ[#, z] & /@ Join[{first}, polys])]];
  ];

  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, y],
    {{Join[{first}, polys]}}
  ]
];

Polyinf[J_, x_, y_] := PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,y],{{{0,0,0,0,0,0,0,0,-(1/21772800)}}}];

LaunchKernels[];

PMP2SDP[datfile_, prec_:600] := Module[
    {
        pols, norm, obj
    },
    pols = 
        Flatten[{
        Flatten[N[ParallelTable[Poly[i, mgap+x, x],{i, 0, 1000, 2}],prec]],
        Flatten[N[ParallelTable[Poly[i, mgap+x, x],{i, 1500, 5000, 100}],prec]],
        Flatten[N[ParallelTable[Poly[i, mgap+x, x],{i, 6000, 20000, 500}],prec]],
        Flatten[N[ParallelTable[Poly[i, mgap+x, x],{i, 20000, 50000, 2000}],prec]],
        Flatten[N[ParallelTable[Poly[i, m1, x],{i, J1, J1, 2}],prec]],
        Flatten[N[ParallelTable[Polyinf[i, mgap+x, x],{i, 0, 0, 2}],prec]]
    },1];

    norm = -1 * N[Flatten[{{0,1},list0}],prec];
	obj  = -1 *N[Flatten[{{1,0},list0}],prec];
    
    Print["size of nomr = ", Length[norm]];
    Print["size of obj = ", Length[obj]];

    WritePmpJson[datfile, SDP[obj, norm, pols], prec, getAnalyticSampleData]
];

PMP2SDP["n_pmp.json", 1000];



