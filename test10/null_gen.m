(* ::Package:: *)

(* null constraint generator for ABAB scattering *)
(* ---------- coefficient helper ---------- *)

g[a_, b_] := gABAB @@ Sort[{a, b}];
c[a_, b_] := cABAB @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* AAAA scattering change the ansatz to be s-u symmetric *)
fwdM[s_, t_, mA_, Nmax_Integer] :=
    gABB^2 * (1/s + 1/t + 1/u) +
    gAAB^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) +
    Total[
        Function[{ab},
            g[ab[[1]], ab[[2]]]
            * (s-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};
    
stuM[s_, t_, mA_, Nmax_Integer] :=
    gABB^2 * (1/s + 1/t + 1/u) +
    gAAB^2 * (1/(s-mA^2) + 1/(t-mA^2) + 1/(u-mA^2)) +
    Total[
        Function[{ab},
            c[ab[[1]], ab[[2]]]
            * (s-2mA^2/3)^ab[[1]]
            * (u-2mA^2/3)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};


(* Test if the double contour still works *)

DblCtrTest[mA_, k_Integer, q_Integer, Nmax_Integer] := Residue[ Residue[ 1/(s*t) ( fwdM[s, t, mA, Nmax]/(s^(k-q)*t^q) - fwdM[t, s, mA, Nmax]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

Print["M(s,t)/((s)^(k-q+1)*(t)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[s,t,mA,10]/(s^(5-2)*t^(2)),{s,Infinity}],{t,0}]//FullSimplify]
Print["M(t,s)/((t)^(k-q+1)*(s)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[t,s,mA,10]/(t^(5-2)*s^(2)),{s,Infinity}],{t,0}]//FullSimplify]


Print["Double contour test: ", DblCtrTest[mA, 5, 2, 10] ];

Print["M(s,t)/((s)^(k-q+1)*(t)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[s,t,mA,10]/(s^(5-3)*t^(3)),{s,Infinity}],{t,0}]//FullSimplify]
Print["M(t,s)/((t)^(k-q+1)*(s)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[t,s,mA,10]/(t^(5-3)*s^(3)),{s,Infinity}],{t,0}]//FullSimplify];
Print["Double contour test: ", DblCtrTest[mA, 5, 3, 10] ];



(* Null constraints *)

(* n4 null constraint from stu low-energy ansatz *)

(* note that c[1,3] =2 c[0,4] and c[2,2] = 3c[0,4] *)
s1[mA_] := 2mA^2/3;
s2[mA_] := 2mA^2-2mA^2/3-t;
Ker[s_,t_,s1_,s2_,k_] := 1/(s-s1)*1/((s-s1)*(s-s2))^(k/2);

Print["g[4, 0] = ", -SeriesCoefficient[Residue[Ker[s,t,s1[m],s2[m],2]*stuM[s,t,m,10],{s,Infinity}],{t,2m^2/3,2}]]

Print["g[4, 0] = ", -SeriesCoefficient[Residue[Ker[s,t,s1[m],s2[m],4]*stuM[s,t,m,10],{s,Infinity}],{t,2m^2/3,0}]]


(* partial-wave decomposition *)
Mhe[s_,t_,m_,J_] := 1/2*s/(s-m^2)*LegendreP[J,1+(2s*t)/(s-m^2)^2];
c41[sp_,m_,J_]:= SeriesCoefficient[(Ker[sp,t,s1[m],s2[m],2]-Ker[2m^2-sp-t,t,s1[m],s2[m],2])*Mhe[sp,t,m,J],{t,2m^2/3,2}];
c42[sp_,m_,J_]:= SeriesCoefficient[(Ker[sp,t,s1[m],s2[m],4]-Ker[2m^2-sp-t,t,s1[m],s2[m],4])*Mhe[sp,t,m,J],{t,2m^2/3,0}];
Print["c[4,0] = ", c41[sp,m,J]//FullSimplify]
Print["c[4,0] = ", c42[sp,m,J]//FullSimplify]

Print["null constraint n4 = ", c41[sp,m,J]-c42[sp,m,J]//FullSimplify]


Print["massless-limit c[4,0] = ", Limit[c41[sp,m,J],{m->0}]//FullSimplify];

Print["massless-limit c[4,0] = ", Limit[c42[sp,m,J],{m->0}]//FullSimplify];

Print["massless-limit null constraint n4 = ", Limit[c41[sp,m,J],{m->0}]-2*Limit[c42[sp,m,J],{m->0}]//FullSimplify];



