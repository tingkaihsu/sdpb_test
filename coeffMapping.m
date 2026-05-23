(* ::Package:: *)

(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

c[a_, b_] := cCoeff @@ Sort[{a, b}];

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
            c[ab[[1]], ab[[2]]]
            * (s-4mA^2/3)^ab[[1]]
            * (u-4mA^2/3)^ab[[2]]
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

fwdKer[sp_, t_, s1_, s2_, k_] := 1/( (sp-s1)*( (sp-s1)*(sp-s2) )^(k/2) );

g2 = -SeriesCoefficient[Residue[fwdKer[sp,t,2m^2,2m^2-t,2]*fwdMlow[sp, t, m, 10], {sp, Infinity}], {t,0,0}]//FullSimplify;

Print["coefficient mapping g2 = ", g2];
Print[""]
g3 = -SeriesCoefficient[Residue[fwdKer[sp,t,2m^2,2m^2-t,2]*fwdMlow[sp,t,m,10], {sp,Infinity}], {t,0,1}]//FullSimplify;

Print["coefficient mapping g3 = ", g3];
Print[""]

c2 = -SeriesCoefficient[Residue[fwdKer[sp,t,2m^2,2m^2-t,2]*stuMlow[sp,t,m,10], {sp, Infinity}], {t,0,0}]//FullSimplify;
Print["coefficient mapping at stu expansion point: c2 = ", c2]



(* Extract the fwd coeff directly *)
Ker[s_, t_, mA_, k_, q_, s1_, s2_] := 1/((s-s1)(4mA^2-s-t-s2)) 1/((s-s1)^(k-q)(4mA^2-s-t-s2)^q);

Print["g[2,0] = ", SeriesCoefficient[Residue[Ker[s,t,mA,2,0,2mA^2,2mA^2]fwdMlow[s,t,mA,10],{s,2mA^2}],{t,0,0}]]
