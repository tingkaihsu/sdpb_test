(* f[sp_, s_, t_, mA_] := ( 1/(sp - s) + 1/(sp - u) )*LegendreP[J, 1+2*t/(sp-4mA^2)]/.{u -> 4mA^2 - s - t}

ma = 0;

g20[x_, J_] := 1/2*Sqrt[ sp/ (sp-4*mA^2) ] * (sp^(-3) + (-4*mA^2 + sp)^(-3))/.{sp -> 1/(1-x), mA -> ma};

g31[x_, J_] := -Sqrt[ sp/ (sp-4*mA^2) ] * ((-3 + J*(1 + J)*(-4*mA^2 + sp)^3*(sp^(-3) + (-4*mA^2 + sp)^(-3)))/(-4*mA^2 + sp)^4)/.{sp -> 1/(1-x), mA -> ma};

Print["g[2, 0] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, mA], {s, 0, 2}], {t, 0, 0}]//FullSimplify];

Print["massless limit g[2, 0] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, 0], {s, 0, 2}], {t, 0, 0}]//FullSimplify];

Print["g[3, 1] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, mA], {s, 0, 2}], {t, 0, 1}]//FullSimplify];

Print["massless limit g[3, 1] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, 0], {s, 0, 2}], {t, 0, 1}]//FullSimplify];

Print["massless limit g[2, 0] = ", g20[x, J]//FullSimplify]

Print["massless limit g[3, 1] = ", g31[x, J]//FullSimplify] *)

(* ============================================================
   AMPLITUDE ANSATZ — corrected
   ============================================================

   BUG FOUND & FIXED:
   ------------------
   validTriples[Nmax] was UNDEFINED in this file, and the intended
   definition used only SORTED pairs {a,b} with a<=b.

   That breaks s<->t symmetry.  For a pair like {0,3}, the sum
       del[0,3] * (s-4/3)^0 * (t-4/3)^3
   only contributes the (t-4/3)^3 term; the mirror (s-4/3)^3
   is never added.  Consequently:

     Mlow[s,t,3] - Mlow[s,t,2]  (buggy)
       = delCoeff[0,3]*(t-4/3)^3 + delCoeff[1,2]*(s-4/3)*(t-4/3)^2
       -- asymmetric; missing the s-channel contributions entirely.

   FIX:
   ----
   validTriples[Nmax] must return ALL ORDERED pairs {a,b} with
   a+b <= Nmax (i.e. both {a,b} and {b,a} for a=/=b).

   Because del[a,b] := delCoeff @@ Sort[{a,b}], the pair {b,a}
   automatically maps to the same symbol delCoeff[min,max] as {a,b}.
   The total contribution for any off-diagonal pair is therefore

       delCoeff[a,b] * ( (s-4/3)^a*(t-4/3)^b
                        +(s-4/3)^b*(t-4/3)^a )

   which is manifestly s<->t symmetric.

   With this fix, Mlow[s,t,3] - Mlow[s,t,2] reduces to PURELY
   degree-3 symmetric terms:

     (-8+3s+3t)/27 * [
       (16 + 9s^2 + 3t(-4+3t) - 3s(4+3t)) * delCoeff[0,3]
      +(-4+3s)(-4+3t)                      * delCoeff[1,2]
     ]
   ============================================================ *)

(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

(* ---------- FIX: validTriples returns ALL ordered pairs ----------
   Returns every {a,b} with a>=0, b>=0, a+b<=Nmax.
   For a=/=b this includes BOTH {a,b} and {b,a}, giving the
   s<->t-symmetric polynomial automatically via the del sort above. *)

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* Quick sanity check (comment out if not needed):
   validTriples[2] should be:
   {{0,0},{0,1},{0,2},{1,0},{1,1},{2,0}}  -- 6 pairs
   validTriples[3] adds:
   {{0,3},{1,2},{2,1},{3,0}}              -- 4 new pairs *)

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference:
   Mlow[s_,t_] := gAAB^2*(1/s+1/t+1/u) + gAAA^2*(1/(s-1)+1/(t-1)+1/(u-1))
                + gAAAA2*(s-4/3)^2 /. {u -> 4-s-t};              *)

ma = 0;

Mlow[s_, t_, mA_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-mA) + 1/(t-mA) + 1/(u-mA)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s)^ab[[1]]
            * (t)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4mA^2 - s - t};

(* Kernel with subtraction points at s1 and s2 *)
KerMlow[sp_, s_, t_, mA_, s1_, s2_, k_, Nmax_Integer] := Mlow[sp, t, mA, Nmax] /(sp) /((sp-s1)*(sp-s2))^(k/2);

s1 = 0;
s2 = -t;

Print["delCoeff[0, 4] - delCoeff[1, 3] + delCoeff[2, 2] = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, s1, s2, 2, 8], {sp, Infinity}], {t, 0, 2}]//FullSimplify];

s1p = 0;
s2p = -t;

Print["delCoeff[0, 4] = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, s1p, s2p, 4, 8], {sp, Infinity}], {t, 0, 0}]//FullSimplify];


g40[sp_, mA_] := SeriesCoefficient[LegendreP[J, 1+2*t/sp]*(1/(sp)/(sp+t)^2/(sp)^2 - 1/(4mA-sp-t)/(4mA-sp-t+t)^2/(4mA-sp-t)^2), {t, 0, 0}];

g22[sp_, mA_] := SeriesCoefficient[LegendreP[J, 1+2*t/sp]*(1/(sp)/(sp+t)/(sp) - 1/(4mA-sp-t)/(4mA-sp-t+t)/(4mA-sp-t)), {t, 0, 2}];

Print["delCoeff[0, 4] = ", g40[sp, mA]//FullSimplify];

Print["delCoeff[0, 4] - delCoeff[1, 3] + delCoeff[2, 2] = ", g22[sp, mA]//FullSimplify];

Print["massless delCoeff[0, 4] = ", g40[sp, 0]//FullSimplify];

Print["massless delCoeff[0, 4] - delCoeff[1, 3] + delCoeff[2, 2] = ", g22[sp, 0]//FullSimplify];

(* differed by a factor of 2, see the following analysis *)

Print["n4 = ", 2*g40[sp, 0] - g22[sp, 0]//FullSimplify];


(* this tells us they should differ by a factor of 2 *)
Print["s^4 coefficient = ", SeriesCoefficient[SeriesCoefficient[(s^2+t^2+u^2)^2/.{u -> -s-t}, {s, 0, 4}], {t, 0, 0}]//FullSimplify ];

Print["s^2 t^2 coefficient = ", SeriesCoefficient[SeriesCoefficient[(s^2+t^2+u^2)^2/.{u -> -s-t}, {s, 0, 2}], {t, 0, 2}]//FullSimplify ];

Print["s t^3 coefficient = ", SeriesCoefficient[SeriesCoefficient[(s^2+t^2+u^2)^2/.{u -> -s-t}, {s, 0, 1}], {t, 0, 3}]//FullSimplify ];


Print["large J limit = ", SeriesCoefficient[Sqrt[sp/(-4*mA^2 + sp)]*(4*mA - sp)^(-5) + (4 - (-2 + J)*J*(1 + J)*(3 + J))/(4*sp^5) + (2*J*(1 + J))/(sp*(-4*mA + sp)^4) - ((-1 + J)*J*(1 + J)*(2 + J))/(4*sp^2*(-4*mA + sp)^3)/.{sp -> 1/(1-x)}, {J, 0, 4}]//FullSimplify];
