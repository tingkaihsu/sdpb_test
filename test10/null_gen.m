(* ::Package:: *)

(* null constraint generator for ABAB scattering *)
(* ---------- coefficient helper ---------- *)
ClearAll["Global`*"];

g[a_, b_] := gABBA @@ Sort[{a, b}];
c[a_, b_] := cAABB @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* ABBA scattering change the ansatz to be s-u symmetric *)
suMABBA[s_, t_, mA_, Nmax_Integer] :=
    gABB^2 * (1/s + 1/u) + gAAB*gBBB * (1/t) + 
    gAAB^2 * (1/(s-mA^2) + 1/(u-mA^2)) + gAAA*gBBA * (1/(t-mA^2)) +
    Total[
        Function[{ab},
            g[ab[[1]], ab[[2]]]
            * (s-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};

(* there is no physical meaning of forward limit for BBAA *)
(* we should use a u <-> t symmetric ansatz*)

tuMAABB[s_, t_, mA_, Nmax_Integer] :=
    gBBB*gAAB * (1/s + 1/u) + gBBA^2 * (1/t) +
    gBBA*gAAA * (1/(s-mA^2) + 1/(u-mA^2)) + gBAA^2 * (1/(t-mA^2)) +
    Total[
        Function[{ab},
            c[ab[[1]], ab[[2]]]
            * (t-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]  
        ] /@ validTriples[Nmax]
    ] /. {u -> 2mA^2 - s - t};

s1[m_,t_] := m^2;
s2[m_,t_] := 4m^2-t-s1[m,t];
Ker[s_, t_, s1_, s2_, k_] := 1/(s-s1)*1/((s-s1)*(s-s2))^(k/2);

Print["g[4,0] = ", SeriesCoefficient[ Residue[Ker[s,t,s1[m,t],s2[m,t],2]*suMABBA[s,t,m,10],{s,s1[m,t]}]+Residue[Ker[s,t,s1[m,t],s2[m,t],2]*suMABBA[s,t,m,10],{s,s2[m,t]}] ,{t,0,2}]//FullSimplify]
Print["g[4,0] = ", SeriesCoefficient[ Residue[Ker[s,t,s1[m,t],s2[m,t],4]*suMABBA[s,t,m,10],{s,s1[m,t]}]+Residue[Ker[s,t,s1[m,t],s2[m,t],4]*suMABBA[s,t,m,10],{s,s2[m,t]}] ,{t,0,0}]//FullSimplify]

(* Note that the 1/2 factor in MABAB *)
MABAB[s_, t_, m_] := 1/2*s/(s-m^2)*LegendreP[J, (-m^4+s (-2 m^2+s+2 t))/(m^2-s)^2];
(* Note that the 1/2 factor in MABBA *)
MABBA[s_,t_,m_] := 1/2*s/(s-m^2)*LegendreP[J, 1+(2 s t)/(m^2-s)^2];
(* BB -> AA scattering is NOT s <-> t symmetric? *)
MBBAA[s_, t_, m_] := (s/(s-4m^2))^(1/4)*LegendreP[J,(-2 m^2+s+2 t)/Sqrt[s (-4 m^2+s)]];

(* ---------- partial-wave kernels around crossing point ---------- *)


(* AB -> AB s-channel *)
sKer[sp_, m_, J_, k_, q_] := Assuming[ J \[Element] Integers && J >= 0 && m >= 0 && sp >= 0, Simplify @ SeriesCoefficient[ 1/((sp-m^2)*t)*1/((sp-m^2)^(k-q)*t^q) * MABBA[sp, t, m] - 1/((sp-m^2)*t)*1/((sp-m^2)^q*t^(k-q)) * MBBAA[sp, t, m], {t, 0, -1} ] ];

(* crossed u-channel contribution *)
uKer[sp_, m_, J_, k_, q_] := Assuming[ J \[Element] Integers && J >= 0 && m >= 0 && sp >= 0, Simplify @ SeriesCoefficient[ 1/(((2m^2-sp-t)-m^2)*t)*1/(((2m^2-sp-t)-m^2)^(k-q)*t^q) * MABBA[sp, t, m] - 1/(((2m^2-sp-t)-m^2)*t)*1/(((2m^2-sp-t)-m^2)^q*t^(k-q)) * MABAB[sp, t, m], {t, 0, -1} ] ];

(* null kernel *)
Xkq[k_Integer, q_Integer, J_, sp_, m_] := Simplify[ sKer[sp, m, J, k, q] - uKer[sp, m, J, k, q] ];
  
  
k1 = 5;
q1 = 2;

Print["X[5,2] = ", Xkq[k1,q1,J,sp,m]//Simplify];

k2 = 4;
q2 = 2;

Print["X[4,2] = ", Xkq[k2,q2,J,sp,m]//Simplify];

(* term by term, n4ABBA *)
n4ABBA[sp_, m_, J_, k_, q_] := Assuming[ J \[Element] Integers && J >= 0 && m >= 0 && sp >= 0, Simplify @ SeriesCoefficient[1/((sp-m^2)*t)*1/((sp-m^2)^(k-q)*t^q) * MABBA[sp, t, m]-1/(((2m^2-sp-t)-m^2)*t)*1/(((2m^2-sp-t)-m^2)^(k-q)*t^q) * MABBA[sp, t, m], {t,0,-1}] ];
n4BBAA[sp_, m_, J_, k_, q_] := Assuming[ J \[Element] Integers && J >= 0 && m >= 0 && sp >= 0, Simplify @ SeriesCoefficient[- (1/((sp-m^2)*t))*1/((sp-m^2)^q*t^(k-q)) * MBBAA[sp, t, m], {t,0,-1}] ];
n4ABAB[sp_, m_, J_, k_, q_] := Assuming[ J \[Element] Integers && J >= 0 && m >= 0 && sp >= 0, Simplify @ SeriesCoefficient[1/(((2m^2-sp-t)-m^2)*t)*1/(((2m^2-sp-t)-m^2)^q*t^(k-q)) * MABAB[sp, t, m], {t,0,-1}] ];

Print["n4ABBA[4,2] = ", n4ABBA[sp,mA,J,k2,q2]//Simplify]
Print["n4BBAA[4,2] = ", n4BBAA[sp,mA,J,k2,q2]//Simplify]
Print["n4ABAB[4,2] = ", n4ABAB[sp,mA,J,k2,q2]//Simplify]



