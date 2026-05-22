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

(* stu symmetric expansion point *)
stuMlow[s_, t_, mA_, Nmax_Integer] :=
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

stuKer[sp_, s_, t_, mA_, s1_, s2_, k_, Nmax_Integer] := stuMlow[sp, t, mA, Nmax] /(sp-4mA^2/3) /((sp-s1)*(sp-s2))^(k/2);

s1 = 4mA^2/3;
s2 = 8mA^2/3-t;

(* --------------------------------------------------------------------------------------------- *)

fwdKer[sp_, s_, t_, mA_, s1_, s2_, k_, Nmax_Integer] := fwdMlow[sp, t, mA, Nmax] / ( (sp-2mA^2)*( (sp-s1)*(sp-s2) )^(k/2) );

(* forward limit *)
s3 = 2mA^2;
s4 = 2mA^2 - t;

(* forward dispersive sum rule *)
fwdKersum[sp_, t_, mA_, s1_, s2_, k_, J_] := ( 1/( (sp-2mA^2)*( (sp-s1)*(sp-s2) )^(k/2) ) - 1/( ( (4mA^2-sp-t)-2mA^2 )*( ( (4mA^2-sp-t)-s1)*((4mA^2-sp-t)-s2) )^(k/2) ) ) * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1+2*t/(sp-4mA^2)];

(* stu dispersive sum rule *)
stuKersum[sp_, t_, mA_, s1_, s2_, k_, J_] := ( 1/( (sp-4mA^2/3)*( (sp-s1)*(sp-s2) )^(k/2) ) - 1/( ( (4mA^2-sp-t) - 4mA^2/3 )*( ((4mA^2-sp-t)-s1)*((4mA^2-sp-t)-s2) )^(k/2) ) ) * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1+2*t/(sp-4mA^2)];

Print[""]

Print["at fwd expansion points and leading order= ", -SeriesCoefficient[Residue[fwdKer[sp, s, t, mA, s3, s4, 2, 8], {sp, Infinity}], {t, 0, 0}]//FullSimplify]

Print["at fwd expansion points and leading order sum rule= ", SeriesCoefficient[fwdKersum[sp, t, mA, s3, s4, 2, J], {t, 0, 0}]//FullSimplify]

G2p = SeriesCoefficient[fwdKersum[sp, t, mA, s3, s4, 2, J], {t, 0, 0}]//FullSimplify;

Print["Above under the massless limit = ", Limit[G2p, mA->0]//FullSimplify]

Print[""]

Print["at fwd expansion points, g[3,1] = ", -SeriesCoefficient[Residue[fwdKer[sp, s, t, mA, s3, s4, 2, 8], {sp, Infinity}], {t, 0, 1}]//FullSimplify]

G3p = SeriesCoefficient[fwdKersum[sp, t, mA, s3, s4, 2, J], {t, 0, 1}]//FullSimplify;

Print["(twice-subtracted) sum rule at fwd = ", SeriesCoefficient[fwdKersum[sp, t, mA, s3, s4, 2, J], {t, 0, 1}]//FullSimplify]

Print["Above under the massless limit = ", Limit[G3p, mA->0]//FullSimplify]

Print[""]