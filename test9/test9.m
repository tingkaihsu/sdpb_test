(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

(* ---------- FIX: validTriples returns ALL ordered pairs ----------
   Returns every {a,b} with a>=0, b>=0, a+b<=Nmax.
   For a=/=b this includes BOTH {a,b} and {b,a}, giving the
   s<->t-symmetric polynomial automatically via the del sort above. *)

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference:
   Mlow[s_,t_] := gAAB^2*(1/s+1/t+1/u) + gAAA^2*(1/(s-1)+1/(t-1)+1/(u-1))
                + gAAAA2*(s-4/3)^2 /. {u -> 4-s-t};              *)

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

(* ============================================================
   Apr 24th Test
   ============================================================ *)

integrand[sp_, t_, s_, s1_, s2_, Nmax_Integer] := Mlow[sp, t, Nmax]/(sp-s)/(sp-s1)/(sp-s2);

(* Apr 25th Test *)
(* k must be even and larger than or equal to 2 *)
Chigh[sp_, t_, s_, s1_, s2_, Nmax_Integer, k_Integer] := Mlow[sp, t, Nmax]/(sp-s)/((sp-s1)*(sp-s2))^(k/2);

C2 = Limit[-Residue[Chigh[sp, t, s, t, 4-t-s, 6, 2], {sp, Infinity}], {s -> 0}];
Print["Nmax = 6, s1 = t, s2 = 4-t-s, k=2: ", C2 // FullSimplify]

C4 = Limit[-Residue[Chigh[sp, t, s, t, 4-s-t, 4, 4], {sp, Infinity}], {s -> 0}];
(* Print["Nmax = 4, s1 = t, s2 = 4-s-t, k=4: ", C4 // FullSimplify] *)

(* UV partial waves *)

f[sp_, s_, t_]:= (1/(sp-t) + 1/(sp-4+2t)) * Sqrt[sp/(sp-4)] * 1/(sp-s)/(sp-u)/.{u->4-s-t};

(* Expansion of Legendre polynomials *)

series[k_Integer] := Series[LegendreP[J, 1+x], {x, 0, k}];

(* Print["LegendreP[J, x] series expansion: ", series[3] // Normal]; *)

g = (1 + (J*(1 + J)*x)/2 + ((-1 + J)*J*(1 + J)*(2 + J)*x^2)/16)*(1/(sp-t) + 1/(sp-4+2t)) * Sqrt[sp/(sp-4)] * 1/(sp-s)/(sp-u)/.{u->4-s-t}/.{x->2t/(sp-4)};

Print["g function: ", SeriesCoefficient[Limit[g, {s->0}], {t, 0, 0}] // FullSimplify];

(* Massless scattering with gravitation have a pole term *)

Mgrav[s_, t_] := 8*Pi*G*( s*t/u + s*u/t + t*u/s ) - d3*( 1/s + 1/t + 1/u ) - d4 + g2*( s^2 + t^2 + u^2 ) + g3 (s*t*u) + g4*( s^2 + t^2 + u^2 )^2/.{u->-s-t};
Cgrav[sp_, t_, s_, s1_, s2_, k_Integer] := Mgrav[sp, t]/(sp-s)/((sp-s1)*(sp-s2))^(k/2);

(* Print["s1 = 0, s2 = -u, k = 2: ", -Residue[Cgrav[sp, t, s, 0, s+t, 2], {sp, Infinity}] // FullSimplify]; *)