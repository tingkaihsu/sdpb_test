(* ::Package:: *)

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

(* Forward amplitude ansatz *)
fwdMlow[s_, t_, mA_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s-2mA^2)^ab[[1]]
            * (t-0)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4mA^2 - s - t};

(* Kernel with subtraction points at s1 and s2 *)
KerMlow[sp_, s_, t_, mA_, s1_, s2_, k_, Nmax_Integer] := Mlow[sp, t, mA, Nmax] /(sp-4mA^2/3) /((sp-s1)*(sp-s2))^(k/2);

s1 = 4mA^2/3;
s2 = 8mA^2/3-t;

Print["(twice-subtracted) stu-symmetric ansatz at s1 = 4mA^2, and s2 = -(t-4mA^2/3) = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, s1, s2, 2, 8], {sp, Infinity}], {t, 4mA^2/3, 2}]//FullSimplify];

s1p = 4mA^2/3;
s2p = 8mA^2/3-t;

Print["(forth-subtracted) stu-symmetric ansatz at s1 = 4mA^2, and s2 = -(t-4mA^2/3) = ", SeriesCoefficient[-Residue[KerMlow[sp, s, t, mA, s1p, s2p, 4, 8], {sp, Infinity}], {t, 4mA^2/3, 0}]//FullSimplify];

(* like the setup (2.3) in Extremal EFT, wich s-4mA^2/3, t-4mA^2/3, and u-4mA^2/3 *)
Mtest[s_, t_, mA_] := gAAB^2 * (1/s + 1/t + 1/u) + gAAA^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) + g0 + 
					  g2*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 ) + g4*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )^2 + 
					  g3*( (s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3) ) + g5*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )*( (s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3) )+
					  g6*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )^3 + g6p*((s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3))^2+g7*( (s-4mA^2/3)^2 + (t-4mA^2/3)^2 + (u-4mA^2/3)^2 )^2*((s-4mA^2/3)*(t-4mA^2/3)*(u-4mA^2/3))/.{u -> 4mA^2-s-t};

KerMtest[sp_, s_, t_, mA_, s1_, s2_, k_] := Mtest[sp, t, mA] /(sp-4mA^2/3) /((sp-s1)*(sp-s2))^(k/2);

s1 = 4mA^2/3;
s2 = 8mA^2/3-t;

Print["(twice-subtracted) test ansatz = ", Series[Limit[-Residue[KerMtest[sp, s, t, mA, s1, s2, 2], {sp, Infinity}],{s->4mA^2/3}], {t, 4mA^2/3, 2}]//FullSimplify//Normal];

Print["(twice-subtracted) test ansatz = ", Series[Limit[-Residue[KerMtest[sp, s, t, mA, s1, s2, 4], {sp, Infinity}],{s->4mA^2/3}], {t, 4mA^2/3, 0}]//FullSimplify//Normal];



