(* ::Package:: *)

(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

(* ---------- FIX: validTriples returns ALL ordered pairs ----------
   Returns every {a,b} with a>=0, b>=0, a+b<=Nmax.
   For a=/=b this includes BOTH {a,b} and {b,a}, giving the
   s<->t-symmetric polynomial automatically via the del sort above. *)

(* This might be a problem since in the forward-limit expansion, we do NOT assume s <-> t symmetry. *)

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

(* stu symmetric expansion point *)
Mlow[s_, t_, mA_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s-4mA^2/3)^ab[[1]]
            * (t-4mA^2/3)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4mA^2 - s - t};

(* Forward amplitude ansatz note that we choose s <-> u expansion points *)
fwdMlow[s_, t_, mA_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (u-2mA^2)^ab[[1]]
            * (s-2mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4mA^2 - s - t};

(* Kernel with subtraction points at s1 and s2 *)

(* KerMlow[sp_, s_, t_, mA_, s1_, s2_, k_, Nmax_Integer] := Mlow[sp, t, mA, Nmax] /(sp-4mA^2/3) /((sp-s1)*(sp-s2))^(k/2);

s1 = 4mA^2/3;
s2 = 8mA^2/3-t;

Print["(twice-subtracted) stu-symmetric ansatz at s1 = 4mA^2, and s2 = -(t-4mA^2/3) = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, s1, s2, 2, 8], {sp, Infinity}], {t, 4mA^2/3, 2}]//FullSimplify];

s1p = 4mA^2/3;
s2p = 8mA^2/3-t;

Print["(fourth-subtracted) stu-symmetric ansatz at s1 = 4mA^2, and s2 = -(t-4mA^2/3) = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, s1p, s2p, 4, 8], {sp, Infinity}], {t, 4mA^2/3, 0}]//FullSimplify]; *)

(* stu sum rule *)

(* KerMsum[sp_, t_, mA_, s1_, s2_, k_, J_] := ( 1/((sp-4mA^2/3)*((sp-s1)*(sp-s2))^(k/2)) - 1/(((4mA^2-sp-t)-4mA^2/3)*(((4mA^2-sp-t)-s1)*((4mA^2-sp-t)-s2))^(k/2)) ) * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1+2*t/(sp-4mA^2)];

Print[""]

Print["at stu expansion points, g[2,0] = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, 4mA^2/3, 8mA^2/3-t, 2, 8], {sp, Infinity}], {t, 4mA^2/3, 0}]//FullSimplify]
Print["(twice-subtracted) sum rule at stu = ", SeriesCoefficient[KerMsum[sp, t, mA, 4mA^2/3, 8mA^2/3-t, 2, J], {t, 4mA^2/3, 0}]//FullSimplify]
G2 = SeriesCoefficient[KerMsum[sp, t, mA, 4mA^2/3, 8mA^2/3-t, 2, J], {t, 4mA^2/3, 0}]//FullSimplify;

Print["Above under the massless limit = ", Limit[G2, mA->0]//FullSimplify]

Print[""]

Print["at stu expansion points, g[3,1] = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, 4mA^2/3, 8mA^2/3-t, 2, 8], {sp, Infinity}], {t, 4mA^2/3, 1}]//FullSimplify]

G3 = SeriesCoefficient[KerMsum[sp, t, mA, 4mA^2/3, 8mA^2/3-t, 2, J], {t, 4mA^2/3, 1}]//FullSimplify;
Print["(twice-subtracted) sum rule at stu = ", SeriesCoefficient[KerMsum[sp, t, mA, 4mA^2/3, 8mA^2/3-t, 2, J], {t, 4mA^2/3, 1}]//FullSimplify]
Print["Above under the massless limit = ", Limit[G3, mA->0]//FullSimplify]

Print[""] *)

(* --------------------------------------------------------------------------------------------- *)

KerfwdMlow[sp_, s_, t_, mA_, s1_, s2_, k_, Nmax_Integer] := fwdMlow[sp, t, mA, Nmax] /(sp-2mA^2) /((sp-s1)*(sp-s2))^(k/2);

(* forward limit *)
s1 = 2mA^2;
s2 = 2mA^2 - t;

Print["(twice-subtracted) forward ansatz at s1 = 2mA^2, and s2 = 2mA^2-t = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, s1, s2, 2, 8], {sp, Infinity}], {t, 0, 2}]//FullSimplify];

Print["(twice-subtracted) forward ansatz at s1 = 2mA^2, and s2 = 2mA^2-t = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, s1, s2, 4, 8], {sp, Infinity}], {t, 0, 0}]//FullSimplify];

(* forward dispersive sum rule *)
KerfwdMsum[sp_, t_, mA_, s1_, s2_, k_, J_] := ( 1/((sp-2mA^2)*((sp-s1)*(sp-s2))^(k/2)) - 1/(((4mA^2-sp-t)-2mA^2)*(((4mA^2-sp-t)-s1)*((4mA^2-sp-t)-s2))^(k/2)) ) * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1+2*t/(sp-4mA^2)];

Print[""]

Print["at fwd expansion points, g[2,0] = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, 2mA^2, 2mA^2-t, 2, 8], {sp, Infinity}], {t, 0, 0}]//FullSimplify]
Print["(twice-subtracted) sum rule at fwd = ", SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 0}]//FullSimplify]

G2p = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 0}]//FullSimplify;

Print["Above under the massless limit = ", Limit[G2p, mA->0]//FullSimplify]

Print[""]

Print["at fwd expansion points, g[3,1] = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, 2mA^2, 2mA^2-t, 2, 8], {sp, Infinity}], {t, 0, 1}]//FullSimplify]
G3p = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 1}]//FullSimplify;
Print["(twice-subtracted) sum rule at fwd = ", SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 1}]//FullSimplify]

Print["Above under the massless limit = ", Limit[G3p, mA->0]//FullSimplify]

Print[""]

Print["at fwd expansion points, g[4,0] = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, 2mA^2, 2mA^2-t, 4, 8], {sp, Infinity}], {t, 0, 0}]//FullSimplify]

G4 = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 4, J], {t, 0, 0}];

Print["at fwd expansion points, the sum rule of g[4,0] = ", G4//FullSimplify];

Print["Above under the massless limit = ", Limit[G4, mA->0]//FullSimplify]

Print[""]

n4 = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 2}] - 2*SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 4, J], {t, 0, 0}]//FullSimplify;

Print["n4 = ", n4]

Print["large J limit of n4 = ", SeriesCoefficient[n4, {J, 0, 4}] // Normal//FullSimplify];

Print[""]

(* like the setup (2.3) in Extremal EFT, wich s-4mA^2/3, t-4mA^2/3, and u-4mA^2/3 *)

(* Mtest[s_, t_, mA_] := gAAB^2 * (1/s + 1/t + 1/u) + gAAA^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) + g0 + 
					  g2*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 ) + g4*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )^2 + 
					  g3*( (s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3) ) + g5*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )*( (s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3) )+
					  g6*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )^3 + g6p*((s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3))^2+g7*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )^2*((s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3))/.{u -> 4mA^2-s-t};

KerMtest[sp_, s_, t_, mA_, s1_, s2_, k_] := Mtest[sp, t, mA] /(sp-4mA^2/3) /((sp-s1)*(sp-s2))^(k/2);

s1 = 4mA^2/3;
s2 = 8mA^2/3-t;

Print["(twice-subtracted) test ansatz = ", Series[Limit[-Residue[KerMtest[sp, s, t, mA, s1, s2, 2], {sp, Infinity}],{s->4mA^2/3}], {t, 4mA^2/3, 1}]//FullSimplify//Normal];

Print["(fourth-subtracted) test ansatz = ", Series[Limit[-Residue[KerMtest[sp, s, t, mA, s1, s2, 4], {sp, Infinity}],{s->4mA^2/3}], {t, 4mA^2/3, 0}]//FullSimplify//Normal]; *)

