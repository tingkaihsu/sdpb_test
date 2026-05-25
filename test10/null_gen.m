(* ::Package:: *)

(* null constraint generator for ABAB scattering *)
(* ---------- coefficient helper ---------- *)

g[a_, b_] := gABAB @@ Sort[{a, b}];
c[a_, b_] := cABAB @@ Sort[{a, b}];
gp[a_, b_] := gBBAA @@ Sort[{a, b}];
cp[a_, b_] := cBBAA @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* ABAB scattering change the ansatz to be s-u symmetric *)
fwdMABAB[s_, t_, mA_, Nmax_Integer] :=
    gABB^2 * (1/s + 1/u) + gAAB*gBBB * (1/t) + 
    gAAB^2 * (1/(s-mA^2) + 1/(u-mA^2)) + gAAA*gBBA * (1/(t-mA^2)) +
    Total[
        Function[{ab},
            g[ab[[1]], ab[[2]]]
            * (s-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};
    
stuMABAB[s_, t_, mA_, Nmax_Integer] :=
    gABB^2 * (1/s + 1/u) + gAAB*gBBB * (1/t) + 
    gAAB^2 * (1/(s-mA^2) + 1/(u-mA^2)) + gAAA*gBBA * (1/(t-mA^2)) +
    Total[
        Function[{ab},
            c[ab[[1]], ab[[2]]]
            * (s-2mA^2/3)^ab[[1]]
            * (u-2mA^2/3)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};

(* there is no physical meaning of forward limit for BBAA *)
fwdMBBAA[s_, t_, mA_, Nmax_Integer] :=
    gBBB*gAAB * (1/s + 1/u) + gBBA^2 * (1/t) +
    gBBA*gAAA * (1/(s-mA^2) + 1/(u-mA^2)) + gBAA^2 * (1/(t-mA^2)) +
    Total[
        Function[{ab},
            gp[ab[[1]], ab[[2]]]
            * (s-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]  
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};

stuMBBAA[s_, t_, mA_, Nmax_Integer] :=
    gBBB*gAAB * (1/s + 1/u) + gBBA^2 * (1/t) +
    gBBA*gAAA * (1/(s-mA^2) + 1/(u-mA^2)) + gBAA^2 * (1/(t-mA^2)) +
    Total[
        Function[{ab},
            cp[ab[[1]], ab[[2]]]
            * (s-2mA^2/3)^ab[[1]]
            * (u-2mA^2/3)^ab[[2]]  
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};