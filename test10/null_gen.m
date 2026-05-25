(* ::Package:: *)

(* null constraint generator for ABAB scattering *)
(* ---------- coefficient helper ---------- *)

g[a_, b_] := gABAB @@ Sort[{a, b}];

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
    ] /. {u -> 4mA^2 - s - t};


(* Test if the double contour still works *)

DblCtrTest[mA_, k_Integer, q_Integer, Nmax_Integer] := Residue[ Residue[ 1/(s*t) ( fwdM[s, t, mA, Nmax]/(s^(k-q)*t^q) - fwdM[t, s, mA, Nmax]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

Print["M(s,t)/((s)^(k-q+1)*(t)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[s,t,mA,10]/(s^(5-2)*t^(2)),{s,Infinity}],{t,0}]//FullSimplify]
Print["M(t,s)/((t)^(k-q+1)*(s)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[t,s,mA,10]/(t^(5-2)*s^(2)),{s,Infinity}],{t,0}]//FullSimplify]


Print["Double contour test: ", DblCtrTest[mA, 5, 2, 10] ];

Print["M(s,t)/((s)^(k-q+1)*(t)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[s,t,mA,10]/(s^(5-3)*t^(3)),{s,Infinity}],{t,0}]//FullSimplify]
Print["M(t,s)/((t)^(k-q+1)*(s)^(q+1)) = ", Residue[Residue[1/(s*t)*fwdM[t,s,mA,10]/(t^(5-3)*s^(3)),{s,Infinity}],{t,0}]//FullSimplify]


Print["Double contour test: ", DblCtrTest[mA, 5, 3, 10] ];

(* Null constraints *)

(* n4 null constraint from stu low-energy ansatz TODO *)

ClearAll[t0,KerstuMsumStable,KerstuCoeffStable,KerstuMsumStableN];

t0[mA_]:=4 mA^2/3;

(*Stable exact rewrite for k=2 or 4*)
KerstuMsumStable[sp_,t_,mA_,k_Integer?((#==2||#==4)&),J_]:=Module[{\[Delta],a,b,c,pref},\[Delta]=t-t0[mA];
a=sp-2 mA^2;
b=sp-t0[mA];
c=sp-2 mA^2/3+\[Delta];
pref=Sqrt[sp/(sp-4 mA^2)]*LegendreP[J,1+2 t/(sp-4 mA^2)];
pref*(1/a+1/c)/(b*(b+\[Delta]))^(k/2)];

(*Series coefficient about t=4/3 mA^2,using the shifted variable \[Delta]*)
KerstuCoeffStable[sp_,mA_,k_Integer?((#==2||#==4)&),J_,n_Integer?NonNegative]:=Module[{\[Delta],expr,assm},expr=KerstuMsumStable[sp,t0[mA]+\[Delta],mA,k,J];
assm=sp>4 mA^2&&0<=mA<=2/5&&Element[J,Integers]&&J>=0;
Assuming[assm,FullSimplify[SeriesCoefficient[expr,{\[Delta],0,n}] ] ] ];

(*High-precision numerical evaluation*)
KerstuMsumStableN[sp_?NumericQ,t_?NumericQ,mA_?NumericQ,k_Integer?((#==2||#==4)&),J_Integer?NonNegative,wp_:80]:=Module[{expr},expr=KerstuMsumStable[SetPrecision[sp,wp],SetPrecision[t,wp],SetPrecision[mA,wp],k,J];
Block[{$MaxExtraPrecision=2 wp},Chop[N[expr,wp] ] ] ];

Print["n4 = ",KerstuCoeffStable[sp,mA,2,J,2]-2*KerstuCoeffStable[sp,mA,4,J,0]//FullSimplify];

(* Large J limit *)
(* Print["Large J limit of n4 = ", Limit[(KerstuCoeffStable[sp,mA,2,J,2]-2*KerstuCoeffStable[sp,mA,4,J,0])/J^4, J -> Infinity]//FullSimplify]; *)

(* take the massless limit *)
Limit[KerstuCoeffStable[sp,mA,2,J,2],mA->0,Direction->"FromAbove",Assumptions->{J\[Element]Integers,sp>0,sp\[Element]Reals}]
Limit[KerstuCoeffStable[sp,mA,4,J,0],mA->0]

(* KerstuCoeffStable[sp,mA,2,J,2]-2*KerstuCoeffStable[sp,mA,4,J,0] *)
(8-8 J-7 J^2+2 J^3+J^4)/(2 sp^5)-2*2/sp^5//FullSimplify



81/4 Sqrt[sp/(-4 mA^2+sp)] (-((8 (22 mA^4-39 mA^2 sp+18 sp^2) LegendreP[J,1+(8 mA^2)/(3 (-4 mA^2+sp))])/((2 mA^2-sp) (8 mA^4-18 mA^2 sp+9 sp^2)^3))+(-2 Sqrt[(8 mA^4-3 mA^2 sp)/(-4 mA^2+sp)^2] (40 mA^4-46 mA^2 sp+9 sp^2) LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+(-8 mA^4+18 mA^2 sp-9 sp^2) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))])/((16 mA^4-14 mA^2 sp+3 sp^2) (8 mA^5-18 mA^3 sp+9 mA sp^2)^2))//FullSimplify
