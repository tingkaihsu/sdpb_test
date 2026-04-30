(* null constraint generator for massless scattering *)


kernelDirect[k_, q_, u_, sp_] :=
  1/u (1/(u^q sp^(k - q + 1)) - 1/(u^(k - q) sp^(q + 1)));

kernelCross[k_, q_, u_, sp_] :=
  1/u (1/(u^q (-sp -u)^(k - q + 1)) - 1/(u^(k - q) (-sp -u)^(q + 1)));



directPiece[k_, q_, J_, t_, sp_] :=
  kernelDirect[k, q, t, u] * LegendreP[J, 1 + 2 t/(sp)]/.{u -> -sp-t};

crossPiece[k_, q_, J_, t_, sp_] :=
  kernelCross[k, q, t, u] * LegendreP[J, 1 + 2 t/(sp)]/.{u -> -sp-t};

directRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[4, 1, J, t, sp], {t, 0, -1}]
];

crossRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[4, 1, J, t, sp], {t, 0, -1}]
];

combined = FullSimplify[directRes - crossRes];

Print["Combined Piece: ", combined];