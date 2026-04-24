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

res3 = Residue[Mlow[s',t, 3]/(s'-s)/(s'-1)/(s'-3+t), {s', s}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-1)/(s'-3+t), {s', 1}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-1)/(s'-3+t), {s', 3-t}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-1)/(s'-3+t), {s', 0}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-1)/(s'-3+t), {s', 4-t}];

Print["res3 = ", res3 // FullSimplify]

res2 = Residue[Mlow[s',t, 2]/(s'-s)/(s'-1)/(s'-3+t), {s', s}] +
      Residue[Mlow[s',t, 2]/(s'-s)/(s'-1)/(s'-3+t), {s', 1}] +
      Residue[Mlow[s',t, 2]/(s'-s)/(s'-1)/(s'-3+t), {s', 3-t}] +
      Residue[Mlow[s',t, 2]/(s'-s)/(s'-1)/(s'-3+t), {s', 0}] +
      Residue[Mlow[s',t, 2]/(s'-s)/(s'-1)/(s'-3+t), {s', 4-t}];

Print["res2 = ", res2 // FullSimplify]

newRes3 = Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', s}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', 1}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', 3-t}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', 0}] +
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', 4-t}] + 
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', t}] + 
      Residue[Mlow[s',t, 3]/(s'-s)/(s'-t)/(s'-4+s+t), {s', 4-s-t}];

Print["newRes3 = ", newRes3 // FullSimplify]

(* UV partial waves *)
f[sp_, s_, t_]:= (1/(sp-s) + 1/(sp-4+s+t)) * Sqrt[sp/(sp-4)] * LegendreP[J, 1 + 2*t/(sp-4)];

