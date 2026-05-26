(* ::Package:: *)

(* null constraint generator for ABAB scattering *)
(* ---------- coefficient helper ---------- *)
ClearAll["Global`*"];

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

Ker[s_, t_, s1_, s2_, k_, q_] := 1/( (s-s1)*(t-s2) ) * 1/(s-s1)^(k-q) * 1/(t-s2)^q;

s1[mA_] := 2mA^2/3;
s2[mA_] := 2mA^2/3;

Print["double contour test (stu ansatz): ", Residue[Residue[ Ker[s,t,s1[m],s2[m],5,2]*stuMABAB[s,t,m,10]-Ker[t,s,s1[m],s2[m],5,2]*stuMBBAA[s,t,m,10], {s,Infinity}], {t,2m^2/3} ]//FullSimplify ]


(* NOT sure about the 1/2 factor in MABAB *)
MABAB[s_, t_, m_] := 1/2*s/(s-m^2)*LegendreP[J, 1+(2*s*t)/(s-m^2)^2];
MBBAA[s_, t_, m_] := (s/(s-m^2))^(1/4)*LegendreP[J,(2t+s-2m^2)/Sqrt[s(s-4m^2)]];

sKer[sp_,m_,J_,k_,q_] := Assuming[J \[Element] Integers && J >= 0,
  Simplify @ SeriesCoefficient[Ker[sp,t,s1[m],s2[m],k,q]*MABAB[sp,t,m]-Ker[t,sp,s1[m],s2[m],k,q]*MBBAA[sp,t,m], {t, 2m^2/3, -1}]
];

uKer[sp_,m_,J_,k_,q_] := Assuming[J \[Element] Integers && J >= 0,
  Simplify @ SeriesCoefficient[Ker[2m^2-sp-t,t,s1[m],s2[m],k,q]*MABAB[sp,t,m]-Ker[t,2m^2-sp-t,s1[m],s2[m],k,q]*MBBAA[sp,t,m], {t, 2m^2/3, -1}]
];

Xkq[k_Integer, q_Integer, J_, sp_, m_] := sKer[sp,m,J,k,q]-uKer[sp,m,J,k,q];

kn = 5;
qn = 2;

Print["X[5,2] = ", Xkq[kn,qn,J,sp,m]//Simplify];
Print["massless-limit X[5,2] = ", Limit[Xkq[kn,qn,J,sp,m],{m -> 0}, Direction->"FromAbove"]];


(4-5*d)J*(J+1)+2(J*(J+1))^2/.{d -> 4}//FullSimplify
(23d^2-12d-20)*J*(J+1)+(-21d-2)*(J*(J+1))^2+4(J*(J+1))^3/.{d->4}//FullSimplify

