ClearAll["Global`*"];

(* String factor without the universal kinematic prefactor *)
StringFactor[s_, t_, u_] := 
  (Gamma[1 - s/4] Gamma[1 - t/4] Gamma[1 - u/4])/
   (Gamma[1 + s/4] Gamma[1 + t/4] Gamma[1 + u/4]);

(* Full type-II four-dilaton amplitude as written in the original file *)
FullAmp[s_, t_] := Module[{u = -s - t},
  (s^2 + t^2 + u^2)^2 * (-1/((s/4) (u/4) (t/4))) * StringFactor[s, t, u]
];

(* Expand uniformly in one low-energy parameter λ.
   This avoids asking for a single monomial coefficient on the slice u = -s - t. *)
LESeries[expr_, order_Integer?NonNegative] := Module[{λ},
  Normal @ Series[expr /. {s -> λ s, t -> λ t}, {λ, 0, order}]
];

(* Low-energy coefficients of the analytic string factor.
   g2 is the λ^2 coefficient (it vanishes).
   g3 is extracted from the λ^3 coefficient, which is proportional to s t u. *)
g2 := FullSimplify[
  SeriesCoefficient[StringFactor[λ s, λ t, -λ (s + t)], {λ, 0, 2}]
];

g3 := FullSimplify[
  SeriesCoefficient[StringFactor[λ s, λ t, -λ (s + t)], {λ, 0, 3}]/
    (s t (-s - t))
];

Print["g2 = ", g2];
Print["g3 = ", g3];

(* If you want the low-energy expansion of the full amplitude itself,
   use:
   LESeries[FullAmp[s, t], n]
   and extract coefficients in λ, not in a single monomial s^n t^0. *)
