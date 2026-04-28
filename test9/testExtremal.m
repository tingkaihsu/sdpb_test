kernelDirect[k_Integer, q_Integer, u_, sp_] :=
  1/u (1/(u^q sp^(k - q + 1)) - 1/(u^(k - q) sp^(q + 1)));

kernelCross[k_Integer, q_Integer, u_, sp_] :=
  1/u (1/(u^q (-sp - u)^(k - q + 1)) - 1/(u^(k - q) (-sp - u)^(q + 1)));



directPiece[k_Integer, q_Integer, J_, t_, sp_] :=
  kernelDirect[k, q, t, sp] * sp * LegendreP[J, 1 + 2 t/sp];

crossPiece[k_Integer, q_Integer, J_, t_, sp_] :=
  kernelCross[k, q, t, sp] * sp * LegendreP[J, 1 + 2 t/sp];

directRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[7, 3, J, t, sp], {t, 0, -1}]
];

crossRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[7, 3, J, t, sp], {t, 0, -1}]
];

combined = FullSimplify[directRes - crossRes];

Print["Combined Piece: ", combined];

(* 5,2 *)
(* -1/36*(J*(1 + J)*(150 + J*(1 + J)*(-43 + 2*J*(1 + J)))) *)

(* 6, 2 *)
(* -1/288*((-3 + J)*J*(1 + J)*(4 + J)*(204 + J*(1 + J)*(-32 + J + J^2))) *)

(* 7, 2 *)
(* -1/14400*(J*(1 + J)*(246960 + J*(1 + J)*(-67908 + J*(1 + J)*(4916 + J*(1 + J)*(-155 + 2*J*(1 + J)))))) *)

(* 7, 3 *)
(* -1/36*((-2 + J)*J*(1 + J)*(3 + J)*(-49 + 2*J*(1 + J))) *)