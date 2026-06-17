(* ::Package:: *)

(* null constraint generator for AAAA scattering *)
(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* AAAA scattering change the ansatz to be s-u symmetric *)
Mfwdlow[s_, t_, mA_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s-2mA^2)^ab[[1]]
            * (u-2mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4mA^2 - s - t};


(* Test if the double contour still works *)

DblCtrTest[mA_, k_Integer, q_Integer, Nmax_Integer] := Residue[ Residue[ 1/(s*t) ( Mfwdlow[s, t, mA, Nmax]/(s^(k-q)*t^q) - Mfwdlow[t, s, mA, Nmax]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

Print["M(s,t)/((s)^(k-q+1)*(t)^(q+1)) = ", Residue[Residue[1/(s*t)*Mfwdlow[s,t,mA,10]/(s^(5-2)*t^(2)),{s,Infinity}],{t,0}]//FullSimplify]
Print["M(t,s)/((t)^(k-q+1)*(s)^(q+1)) = ", Residue[Residue[1/(s*t)*Mfwdlow[t,s,mA,10]/(t^(5-2)*s^(2)),{s,Infinity}],{t,0}]//FullSimplify]


Print["Double contour test: ", DblCtrTest[mA, 5, 2, 10] ];

Print["M(s,t)/((s)^(k-q+1)*(t)^(q+1)) = ", Residue[Residue[1/(s*t)*Mfwdlow[s,t,mA,10]/(s^(5-3)*t^(3)),{s,Infinity}],{t,0}]//FullSimplify]
Print["M(t,s)/((t)^(k-q+1)*(s)^(q+1)) = ", Residue[Residue[1/(s*t)*Mfwdlow[t,s,mA,10]/(t^(5-3)*s^(3)),{s,Infinity}],{t,0}]//FullSimplify]


Print["Double contour test: ", DblCtrTest[mA, 5, 3, 10] ];

(* Null constraints *)

kernelDirect[k_Integer, q_Integer, u_, sp_, mA_] :=
  1/u (1/(u^q sp^(k - q + 1)) - 1/(u^(k - q) sp^(q + 1)));

kernelCross[k_Integer, q_Integer, u_, sp_, mA_] :=
  1/u (1/(u^q (4mA^2 -sp -u)^(k - q + 1)) - 1/(u^(k - q) (4mA^2 -sp -u)^(q + 1)));



directPiece[k_Integer, q_Integer, J_, t_, sp_, mA_] :=
  kernelDirect[k, q, t, sp, mA] * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1 + 2 t/(sp-4mA^2)];

crossPiece[k_Integer, q_Integer, J_, t_, sp_, mA_] :=
  kernelCross[k, q, t, sp, mA] * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1 + 2 t/(sp-4mA^2)];

directRes[k_Integer, q_Integer] := Assuming[J \[Element] Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[k, q, J, t, sp, mA], {t, 0, -1}]
];

crossRes[k_Integer, q_Integer] := Assuming[J \[Element] Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[k, q, J, t, sp, mA], {t, 0, -1}]
];


Print["X52 differs X53 by a minus sign..."]
kn = 5;
qn = 2;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece: ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^6, {J -> Infinity}]//FullSimplify];

kn = 5;
qn = 3;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece: ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^6, {J -> Infinity}]//FullSimplify];


Print["X62 differs X63 by a minus sign..."]
kn = 6;
qn = 2;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece: ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^8, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 6;
qn = 3;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece: ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^8, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 6;
qn = 4;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece: ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^8, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];


kn = 7;
qn = 2;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^10, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 7;
qn = 3;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^10, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 7;
qn = 4;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^10, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 7;
qn = 5;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^10, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];


kn = 8;
qn = 2;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^12, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 8;
qn = 3;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^12, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 8;
qn = 4;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^12, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 8;
qn = 5;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^12, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];

kn = 8;
qn = 6;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q] ];

Print["Combined Piece X", kn, qn, " : ", combined[kn, qn] ];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^12, {J -> Infinity},Assumptions->{sp>0,mA>0}]//FullSimplify];


Limit[ ((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (5+J) (7+J))+((-1+J) (2+J) (-4 mA^2+sp)^4 (64 mA^2+(-28+J+J^2) sp))/sp^5))/(576 (-4 mA^2+sp)^8),mA->0]//FullSimplify

Limit[  -((Sqrt[sp/(-4 mA^2+sp)] (-32 mA^6+24 mA^4 sp-6 mA^2 sp^2+sp^3))/(288 sp^3 (-4 mA^2+sp)^7)),mA->0]//FullSimplify



