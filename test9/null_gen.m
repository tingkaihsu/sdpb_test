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

Print["Double contour test: ", DblCtrTest[mA, 5, 2, 10]];

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

kn = 5;
qn = 2;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q]];

Print["Combined Piece: ", combined[kn, qn]];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^6, {J -> Infinity}]//FullSimplify];

kn = 5;
qn = 3;
combined[k_Integer, q_Integer] := FullSimplify[directRes[k, q] - crossRes[k, q]];

Print["Combined Piece: ", combined[kn, qn]];
(* large J limit *)

Print["Large J limit: ", Limit[combined[kn, qn]/J^6, {J -> Infinity}]//FullSimplify];