(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

Mlow[s_, t_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-1) + 1/(t-1) + 1/(u-1)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s - 4/3)^ab[[1]]
            * (t - 4/3)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4 - s - t};

integrand[sp_, t_, s_, s1_, s2_, Nmax_Integer] := Mlow[sp, t, Nmax]/(sp-s)/(sp-s1)/(sp-s2);

(* Apr 25th Test *)
(* k must be even and larger than or equal to 2 *)
Chigh[sp_, t_, s_, s1_, s2_, Nmax_Integer, k_Integer] := Mlow[sp, t, Nmax]/(sp-s)/((sp-s1)*(sp-s2))^(k/2);

C2 = Limit[-Residue[Chigh[sp, t, s, t, 4-s-t, 6, 2], {sp, Infinity}], {s -> 0}];
Print["Nmax = 6, s1 = t, s2 = 4-s-t, k=2: ", C2 // FullSimplify]

(* result: delCoeff[0, 2] + (16/3 + (-4 + t)*t)*delCoeff[0, 4] + ((-4 + 3*t)*(3*delCoeff[1, 2] + (-4 + 3*t)*delCoeff[2, 2]))/9 *)
(* 16/9 - 4/3 * (t-4/3) +  (t-4/3)^2 *)

C4 = Limit[-Residue[Chigh[sp, t, s, t, 4-s-t, 6, 4], {sp, Infinity}], {s -> 0}];
(* Print["Nmax = 6, s1 = t, s2 = 4-s-t, k=4: ", C4 // FullSimplify] *)

Print["Improved, Nmax = 6, s1 = t, s2 = 4-s-t, k=2: ", C2 - ( 16/9 - 4/3*(t-4/3) )*C4 // FullSimplify]

(* UV partial waves *)

f[sp_, s_, t_]:= (1/(sp-t) + 1/(sp-4+2t)) * Sqrt[sp/(sp-4)] * 1/(sp-s)/(sp-u)/.{u->4-s-t};

(* Expansion of Legendre polynomials *)

series[k_Integer] := Series[LegendreP[J, 1+x], {x, 0, k}];

(* Print["LegendreP[J, x] series expansion: ", series[3] // Normal]; *)

g = (1 + (J*(1 + J)*x)/2 + ((-1 + J)*J*(1 + J)*(2 + J)*x^2)/16)*(1/(sp-t) + 1/(sp-4+2t)) * Sqrt[sp/(sp-4)] * 1/(sp-s)/(sp-u)/.{u->4-s-t}/.{x->2t/(sp-4)};

(* Print["g function: ", SeriesCoefficient[Limit[g, {s->0}], {t, 0, 0}] // FullSimplify]; *)