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


(* dispersion representation of Wilson coefficients *)

v[l_, q_] := Product[l*(l + 1) - a*(a - 1), {a, 1, q}] / (Factorial[q])^2;

Print["v[l,0] = ", v[l, 0]//FullSimplify];
(* 1 *)

Print["v[l,1] = ", v[l, 1]//FullSimplify];
(* l(l+1) *)

(* let the mass be m = 0.2 so that 4m^2 < M^2 = 1 where M  = 1 to infinity *)

ma = 0.2;

f1[x_?NumericQ, J_?IntegerQ] := Sqrt[ M^2/ (M^2-4*mA^2) ] * v[J, 0] * (1/M^2)^2/.{M^2 -> 1/(1-x), mA -> ma};

f1[x_?NumericQ, J_?IntegerQ] := Sqrt[ M^2/ (M^2-4*mA^2) ] * v[J, 1] * (1/M^2)^2 * 1/(M^2-4*mA^2)/.{M^2 -> 1/(1-x), mA -> ma};

(* Null constraints *)

Jmax = 40;

null1[q_Integer, t_, J, sp_, mA_] := 1/t^(q+1) * LegendreP[J, 1 + 2 * t / (sp - 4mA^2)];

(* Print["null1[0] = ", Normal @ Series[null1[0, t, J, sp, mA], {t, 0, 0}]//FullSimplify]; *)

Print["null1[0] = ", Assuming[J ∈ Integers && J >= 0,
  FullSimplify[
    SeriesCoefficient[null1[0, t, J, sp, mA], {t, 0, -1}]
  ]
]];
(* 1 *)

null2[k_Integer, q_Integer, t_, J, sp_, mA_] := -1/t^(k-q+1) * LegendreP[J, 1 + 2 * t / (sp - 4mA^2)];

Print["null2[4, 0] = ", Assuming[J ∈ Integers && J >= 0, Normal @ SeriesCoefficient[null2[4, 0, t, J, sp, mA], {t, 0, -1}]]//FullSimplify];
(* -1/576*((-3 + J)*(-2 + J)*(-1 + J)*J*(1 + J)*(2 + J)*(3 + J)*(4 + J))/(-4*mA^2 + sp)^4 *)

(* u channel cut *)
null3[k_Integer, q_Integer, t_, J, sp_, mA_] := 1/t^(q+1) * 1/(4mA^2-sp-t)^(k-q+1) * LegendreP[J, 1 + 2 * t / (sp - 4mA^2)];

Print["null3[4, 0] = ", Assuming[J ∈ Integers && J >= 0, Normal @ SeriesCoefficient[null3[4, 0, t, J, sp, mA], {t, 0, -1}]]//FullSimplify];
(* (4*mA^2 - sp)^(-5) *)

null4[k_Integer, q_Integer, t_, J, sp_, mA_] := -1/t^(k-q+1) * 1/(4mA^2-sp-t)^(q+1) * LegendreP[J, 1 + 2 * t / (sp - 4mA^2)];

Print["null4[4, 0] = ", Assuming[J ∈ Integers && J >= 0, Normal @ SeriesCoefficient[null4[4, 0, t, J, sp, mA], {t, 0, -1}]]//FullSimplify];
(* -1/576*(576 + (-4 + J)*(-2 + J)*J*(1 + J)*(3 + J)*(5 + J)*(-10 + J + J^2))/(4*mA^2 - sp)^5 *)

a = 1;

b[J_, mA_, sp_] := -1/576*((-3 + J)*(-2 + J)*(-1 + J)*J*(1 + J)*(2 + J)*(3 + J)*(4 + J))/(-4*mA^2 + sp)^4;

c[J_, mA_, sp_] := (4*mA^2 - sp)^(-5);

d[J_, mA_, sp_] := -1/576*(576 + (-4 + J)*(-2 + J)*J*(1 + J)*(3 + J)*(5 + J)*(-10 + J + J^2))/(4*mA^2 - sp)^5;

null[k_Integer, q_Integer, J_, mA_, sp_] = a / (sp)^(k-q+1) + b[J, mA, sp] / (sp)^(q+1) - c[J, mA, sp] + d[J, mA, sp];

Print["null[4, 0] = ", null[4, 0, J, mA, sp]//FullSimplify];

(* null constraint for k = 4, q = 0 *)
(* sp^(-5) + 2/(-4*mA^2 + sp)^5 - ((-2 + J)*J*(1 + J)*(3 + J)*((-3 + J)*(-1 + J)*(2 + J)*(4 + J)*mA^2 - 4*(-11 + J + J^2)*sp))/(144*(4*mA^2 - sp)^5*sp) *)
