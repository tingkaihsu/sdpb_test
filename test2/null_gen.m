(* null constraint generator for massless scattering *)


kernelDirect[k_Integer, q_Integer, u_, sp_] :=
  1/u (1/(u^q sp^(k - q + 1)) - 1/(u^(k - q) sp^(q + 1)));

kernelCross[k_Integer, q_Integer, u_, sp_] :=
  1/u (1/(u^q (-sp -u)^(k - q + 1)) - 1/(u^(k - q) (-sp -u)^(q + 1)));



directPiece[k_Integer, q_Integer, J_, t_, sp_] :=
  kernelDirect[k, q, t, sp] * LegendreP[J, 1 + 2 t/(sp)];

crossPiece[k_Integer, q_Integer, J_, t_, sp_] :=
  kernelCross[k, q, t, sp] * LegendreP[J, 1 + 2 t/(sp)];

directRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[9, 2, J, t, sp], {t, 0, -1}]
];

crossRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[9, 2, J, t, sp], {t, 0, -1}]
];

combined = FullSimplify[directRes - crossRes];

Print["Combined Piece: ", combined];

(* to be safe, we should ONLY use q >= 2, and k - q >= 2 -> k >= 4 *)
(* k, q *)
(* 5, 2 *)
(* -1/36*(J*(1 + J)*(150 + J*(1 + J)*(-43 + 2*J*(1 + J)))) *)

(* 6, 2 *)
(* -1/288*((-3 + J)*J*(1 + J)*(4 + J)*(204 + J*(1 + J)*(-32 + J + J^2))) *)

(* 7, 2 *)
(* -1/14400*(J*(1 + J)*(246960 + J*(1 + J)*(-67908 + J*(1 + J)*(4916 + J*(1 + J)*(-155 + 2*J*(1 + J)))))) *)

(* 8, 2 *)
(* -1/259200*(J*(1 + J)*(-6808320 + J*(1 + J)*(1906416 + J*(1 + J)*(-170976 + J*(1 + J)*(6568 + J*(1 + J)*(-124 + J + J^2)))))) *)

(* 9, 2 *)
(* -1/25401600*(J*(1 + J)*(1015701120 + J*(1 + J)*(-306848736 + J*(1 + J)*(28977336 + J*(1 + J)*(-1293996 + J*(1 + J)*(30170 + J*(1 + J)*(-371 + 2*J*(1 + J)))))))) *)

(* large J limit *)
X52[x_, J_] := -1/36*(J*(1 + J)*(150 + J*(1 + J)*(-43 + 2*J*(1 + J))));
X62[x_, J_] := -1/288*((-3 + J)*J*(1 + J)*(4 + J)*(204 + J*(1 + J)*(-32 + J + J^2)));
X72[x_, J_] := -1/14400*(J*(1 + J)*(246960 + J*(1 + J)*(-67908 + J*(1 + J)*(4916 + J*(1 + J)*(-155 + 2*J*(1 + J))))));
X82[x_, J_] := -1/259200*(J*(1 + J)*(-6808320 + J*(1 + J)*(1906416 + J*(1 + J)*(-170976 + J*(1 + J)*(6568 + J*(1 + J)*(-124 + J + J^2))))));
X92[x_, J_] := -1/25401600*(J*(1 + J)*(1015701120 + J*(1 + J)*(-306848736 + J*(1 + J)*(28977336 + J*(1 + J)*(-1293996 + J*(1 + J)*(30170 + J*(1 + J)*(-371 + 2*J*(1 +	J))))))));

Print["Large J limit: ", SeriesCoefficient[X52[x, J], {J, 0, 14}]//FullSimplify];