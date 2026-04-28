(* null constraint generator for AAAA scattering *)
(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* AAAA scattering *)
Mlow[s_, t_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-1) + 1/(t-1) + 1/(u-1)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s)^ab[[1]]
            * (t)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4 - s - t};


(* Print["v[l,0] = ", v[l, 0]//FullSimplify]; *)
(* 1 *)

(* Print["v[l,1] = ", v[l, 1]//FullSimplify]; *)
(* l(l+1) *)

(* Null constraints *)

kernelDirect[k_Integer, q_Integer, u_, sp_, mA_] :=
  1/u (1/(u^q sp^(k - q + 1)) - 1/(u^(k - q) sp^(q + 1)));

kernelCross[k_Integer, q_Integer, u_, sp_, mA_] :=
  1/u (1/(u^q (4mA^2 -sp -u)^(k - q + 1)) - 1/(u^(k - q) (4mA^2 -sp -u)^(q + 1)));



directPiece[k_Integer, q_Integer, J_, t_, sp_, mA_] :=
  kernelDirect[k, q, t, sp, mA] * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1 + 2 t/(sp-4mA^2)];

crossPiece[k_Integer, q_Integer, J_, t_, sp_, mA_] :=
  kernelCross[k, q, t, sp, mA] * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1 + 2 t/(sp-4mA^2)];

directRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[5, 2, J, t, sp, mA], {t, 0, -1}]
];

crossRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[5, 2, J, t, sp, mA], {t, 0, -1}]
];

combined = FullSimplify[directRes - crossRes];

Print["Combined Piece: ", combined];

(* k, q *)
(* 5, 2 *)
(* (J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6) *)

(* 6, 2 *)
(* (J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(-2304*mA^4 + 1152*mA^2*sp + (-72 + J*(1 + J)*(-18 + J + J^2))*sp^2))/sp^5))/(576*(-4*mA^2 + sp)^7) *)

(* 7, 2 *)
(* (J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(-47 + J + J^2))/(-4*mA^2 + sp)^8) + ((-1 + J)*(2 + J)*(((-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*sp^3)/(4*mA^2 - sp)^5 + 3600/(-4*mA^2 + sp)^2))/sp^6))/14400 *)

(* 8, 2 *)
(* (J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(((-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(-38 + J + J^2))/(4*mA^2 - sp)^9 + ((-1 + J)*(2 + J)*(-(((-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*sp^4)/(-4*mA^2 + sp)^6) + 129600/(-4*mA^2 + sp)^2))/sp^7))/518400 *)
