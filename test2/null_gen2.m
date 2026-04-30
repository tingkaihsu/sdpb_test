(* null constraint generator for massless scattering *)

(* t -> u in the second term *)
kernelDirect[k_, q_, t_, sp_] :=
  1/(t^(q + 1) sp^(k - q + 1)) - 1/((sp)^(k - q + 1) (-sp-t)^(q + 1));

(* u channel cut *)
kernelCross[k_, q_, t_, sp_] :=
  1/(t^(q + 1) (-sp-t)^(k - q + 1)) - 1/((-sp-t)^(k - q + 1) (sp)^(q + 1));



directPiece[k_, q_, J_, t_, sp_] :=
  kernelDirect[k, q, t, sp] * LegendreP[J, 1 + 2 t/(sp)];

crossPiece[k_, q_, J_, t_, sp_] :=
  kernelCross[k, q, t, sp] * LegendreP[J, 1 + 2 t/(sp)];

directRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[5, 2, J, t, sp], {t, 0, -1}]
];

crossRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[5, 2, J, t, sp], {t, 0, -1}]
];

combined = FullSimplify[directRes - crossRes];

Print["Combined Piece: ", combined];