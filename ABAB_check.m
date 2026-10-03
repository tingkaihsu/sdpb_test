(* ::Package:: *)

TestABAB[s_, t_, m_] := gab^2/(s-m^2)+gab^2/(u-m^2)+(gaa gbb)/(t-m^2)/.{u->2m^2-s-t};

TestAABB[s_, t_, m_] := (gaa gbb)/(s-m^2)+gab^2/(u-m^2)+gab^2/(t-m^2)/.{u->2m^2-s-t};

TestAABB[t,s,m]-TestABAB[s,t,m]


gABABCoeff[a_, b_] := gABAB @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];
    
MABAB[s_, t_, mA_, Nmax_Integer] :=
    Total[
        Function[{ab},
            gABABCoeff[ab[[1]], ab[[2]]]
            * (s-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 2*mA^2 - s - t};
    

gAABBCoeff[a_,b_] := gABAB @@ Sort[{a,b}];

MAABB[s_,t_,m_,Nmax_Integer] :=
	Total[
		Function[{ab},
			gAABBCoeff[ab[[1]],ab[[2]]]
			* (t-m^2)^ab[[1]]
			* (u-m^2)^ab[[2]]
			] /@ validTriples[Nmax]
		]/.{u->2m^2-s-t};
		
MABAB[s,t,m,10]-MAABB[t,s,m,10]

(* test Kernel *)

KerLow[k_,q_,Nmax_Integer]:= SeriesCoefficient[SeriesCoefficient[MABAB[s,t,m,Nmax]/((s-m^2)^(k-q+1)t^(q+1)),{s,m^2,-1}],{t,0,-1}]-SeriesCoefficient[SeriesCoefficient[MAABB[s,t,m,Nmax]/((t-m^2)^(k-q+1)s^(q+1)),{s,0,-1}],{t,m^2,-1}];
KerLow[2,1,10]


(* sum rule *)

(* partial waves *)
PartialWaveD[d_, J_, z_] := Hypergeometric2F1[-J, J + d - 3, (d - 2)/2, (1 - z)/2];
phiABAB[s_,m_]:= (s - m^2)^7/s^4;
phiAABB[s_,m_]:= s^(5/4) (s - 4 m^2)^(7/4);

KerABAB[x_,J_,k_,q_,m_]:=Residue[(1/((x-m^2)^(k-q+1)*t^(q+1))-1/((2m^2-x-t-m^2)^(k-q+1)*t^(q+1)))*phiABAB[x,m]*PartialWaveD[10,J,1 + 2 x t/(x - m^2)^2],{t,0}];

(* term 2: s' on the neutral cut ... *)
KerAABB1[x_,J_,k_,q_,m_]:=Residue[(1/((x)^(q+1)*(t-m^2)^(k-q+1)))*phiAABB[x,m]*PartialWaveD[10,J,(2 t + x - 2 m^2)/Sqrt[x (x - 4 m^2)]],{t,m^2}];

(* term 3: the second term on the u-channel cut *)
KerAABB2[x_,J_,k_,q_,m_]:=Residue[-(1/((2m^2-x-t)^(q+1)*(t-m^2)^(k-q+1)))*phiABAB[x,m]*(-1)^J*PartialWaveD[10,J,-1-(2 (2 m^2-t-x) x)/(-m^2+x)^2],{t,m^2}];


NABAB[x_,m_,J_,k_, q_] := KerABAB[x, J, k, q, m] - KerAABB2[x, J, k, q, m];
NAABB[x_,m_,J_,k_, q_] := KerAABB1[x, J, k, q, m];


abab00 = NABAB[x,m,J,0,0]//FullSimplify;
aabb00 = NAABB[x,m,J,0,0]//FullSimplify;

abab10 = NABAB[x,m,J,1,0]//FullSimplify;
aabb10 = NAABB[x,m,J,1,0]//FullSimplify;

abab11 = NABAB[x,m,J,1,1]//FullSimplify;
aabb11 = NAABB[x,m,J,1,1]//FullSimplify;


abab20 = NABAB[x,m,J,2,0]//FullSimplify;
aabb20 = NAABB[x,m,J,2,0]//FullSimplify;

abab21 = NABAB[x,m,J,2,1]//FullSimplify;
aabb21 = NAABB[x,m,J,2,1]//FullSimplify;

abab22 = NABAB[x,m,J,2,2]//FullSimplify;
aabb22 = NAABB[x,m,J,2,2]//FullSimplify;


abab30 = NABAB[x,m,J,3,0]//FullSimplify;
aabb30 = NAABB[x,m,J,3,0]//FullSimplify;

abab31 = NABAB[x,m,J,3,1]//FullSimplify;
aabb31 = NAABB[x,m,J,3,1]//FullSimplify;

abab32 = NABAB[x,m,J,3,2]//FullSimplify;
aabb32 = NAABB[x,m,J,3,2]//FullSimplify;

abab33 = NABAB[x,m,J,3,3]//FullSimplify;
aabb33 = NAABB[x,m,J,3,3]//FullSimplify;


abab40 = NABAB[x,m,J,4,0]//FullSimplify;
aabb40 = NAABB[x,m,J,4,0]//FullSimplify;

abab41 = NABAB[x,m,J,4,1]//FullSimplify;
aabb41 = NAABB[x,m,J,4,1]//FullSimplify;

abab42 = NABAB[x,m,J,4,2]//FullSimplify;
aabb42 = NAABB[x,m,J,4,2]//FullSimplify;

abab43 = NABAB[x,m,J,4,3]//FullSimplify;
aabb43 = NAABB[x,m,J,4,3]//FullSimplify;

abab44 = NABAB[x,m,J,4,4]//FullSimplify;
aabb44 = NAABB[x,m,J,4,4]//FullSimplify;

abab50 = NABAB[x,m,J,5,0]//FullSimplify;
aabb50 = NAABB[x,m,J,5,0]//FullSimplify;

abab51 = NABAB[x,m,J,5,1]//FullSimplify;
aabb51 = NAABB[x,m,J,5,1]//FullSimplify;

abab52 = NABAB[x,m,J,5,2]//FullSimplify;
aabb52 = NAABB[x,m,J,5,2]//FullSimplify;

abab53 = NABAB[x,m,J,5,3]//FullSimplify;
aabb53 = NAABB[x,m,J,5,3]//FullSimplify;

abab54 = NABAB[x,m,J,5,4]//FullSimplify;
aabb54 = NAABB[x,m,J,5,4]//FullSimplify;

abab55 = NABAB[x,m,J,5,5]//FullSimplify;
aabb55 = NAABB[x,m,J,5,5]//FullSimplify;

abab60 = NABAB[x,m,J,6,0]//FullSimplify;
aabb60 = NAABB[x,m,J,6,0]//FullSimplify;

abab61 = NABAB[x,m,J,6,1]//FullSimplify;
aabb61 = NAABB[x,m,J,6,1]//FullSimplify;

abab62 = NABAB[x,m,J,6,2]//FullSimplify;
aabb62 = NAABB[x,m,J,6,2]//FullSimplify;

abab63 = NABAB[x,m,J,6,3]//FullSimplify;
aabb63 = NAABB[x,m,J,6,3]//FullSimplify;

abab64 = NABAB[x,m,J,6,4]//FullSimplify;
aabb64 = NAABB[x,m,J,6,4]//FullSimplify;

abab65 = NABAB[x,m,J,6,5]//FullSimplify;
aabb65 = NAABB[x,m,J,6,5]//FullSimplify;

abab66 = NABAB[x,m,J,6,6]//FullSimplify;
aabb66 = NAABB[x,m,J,6,6]//FullSimplify;

abab70 = NABAB[x,m,J,7,0]//FullSimplify;
aabb70 = NAABB[x,m,J,7,0]//FullSimplify;

abab71 = NABAB[x,m,J,7,1]//FullSimplify;
aabb71 = NAABB[x,m,J,7,1]//FullSimplify;

abab72 = NABAB[x,m,J,7,2]//FullSimplify;
aabb72 = NAABB[x,m,J,7,2]//FullSimplify;

abab73 = NABAB[x,m,J,7,3]//FullSimplify;
aabb73 = NAABB[x,m,J,7,3]//FullSimplify;

abab74 = NABAB[x,m,J,7,4]//FullSimplify;
aabb74 = NAABB[x,m,J,7,4]//FullSimplify;

abab75 = NABAB[x,m,J,7,5]//FullSimplify;
aabb75 = NAABB[x,m,J,7,5]//FullSimplify;

abab76 = NABAB[x,m,J,7,6]//FullSimplify;
aabb76 = NAABB[x,m,J,7,6]//FullSimplify;

abab77 = NABAB[x,m,J,7,7]//FullSimplify;
aabb77 = NAABB[x,m,J,7,7]//FullSimplify;

abab74


lstabab = {abab00, abab10, abab11, abab20, abab21, abab22, abab30, abab31, abab32, abab33, abab40, abab41, abab42, abab43, abab44, abab50, abab51, abab52, abab53, abab54, abab55, abab60, abab61, abab62, abab63, abab64, abab65, abab66, abab70, abab71, abab72, abab73, abab74, abab75, abab76, abab77};

(* safe copy *)
ToString[lstabab, InputForm, PageWidth->Infinity]


lstaabb = {aabb00, aabb10, aabb11, aabb20, aabb21, aabb22, aabb30, aabb31, aabb32, aabb33, aabb40, aabb41, aabb42, aabb43, aabb44, aabb50, aabb51, aabb52, aabb53, aabb54, aabb55, aabb60, aabb61, aabb62, aabb63, aabb64, aabb65, aabb66, aabb70, aabb71, aabb72, aabb73, aabb74, aabb75, aabb76, aabb77}


Length[lstabab]

