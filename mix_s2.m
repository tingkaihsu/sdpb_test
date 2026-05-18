(* ::Package:: *)

(* Exchanging spin-2 massive heavy state *)

Print["LegendreP[2, 1+\!\(\*FractionBox[\(2  t\), \(s - 4  m^2\)]\)]= ", LegendreP[2,1+(2t)/(s-4m^2)]//FullSimplify];

Ms2[s_, u_, M_, m_, g_:1] := g^2/M^4*((t - u)^2 - (4 m^2 - s)^2/3)/(s-M^2)+g^2/M^4*((t - s)^2 - (4 m^2 - u)^2/3)/(u-M^2)+g^2/M^4*((u - s)^2 - (4 m^2 - t)^2/3)/(t-M^2)/.{t -> 4m^2-s-u};


Print["\!\(\*FractionBox[\(s - channel\\\ amplitude\), \(LegendreP[2, \\\ 1 + \*FractionBox[\(2  t\), \(s - 4  m^2\)]]\)]\)= ", (g^2/M^4 ((t - u)^2 - (4 m^2 - s)^2/3)/(s-M^2))/LegendreP[2,1+(2t)/(s-4m^2)]/.{u -> 4m^2-s-t}//FullSimplify]

Print["s-u symmetric = ", Ms2[s,u,M,m,1]-Ms2[u,s,M,m,1]//FullSimplify]


(* coefficient *)
Print["Low-energy ansatz = g[k, q](s-2m^2)^(k-q)(u-2m^2)^q"]

Print["g[3,0]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1], {s, 2m^2, 3}], {u, 2m^2, 0}]//FullSimplify ];
Print["g[3,3]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1], {s, 2m^2, 0}], {u, 2m^2, 3}]//FullSimplify ];

Print["g[3,2]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1], {s, 2m^2, 1}], {u, 2m^2, 2}]//FullSimplify ];
Print["g[3,1]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1], {s, 2m^2, 2}], {u, 2m^2, 1}]//FullSimplify ];

Print["g[2,0]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1],{s, 2m^2, 2}],{u, 2m^2, 0}]//FullSimplify ];
Print["g[2,1]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1],{s, 2m^2, 1}],{u, 2m^2, 1}]//FullSimplify ];
Print["g[2,2]= ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1],{s, 2m^2, 0}],{u, 2m^2, 2}]//FullSimplify ];


 g33[m_,M_]:=  (2 (-8 m^4+4 m^2 M^2+M^4-(M^8 (-8 m^4+4 m^2 M^2+M^4))/(-2 m^2+M^2)^4))/(3 M^12);
 
 g32[m_,M_]:=(2 (-8 m^4+4 m^2 M^2+M^4 (-1+(4 (-m^2 M^4+M^6))/(2 m^2-M^2)^3)))/M^12;
 g22[m_,M_]:= -((4 (32 m^10-64 m^8 M^2+44 m^6 M^4-2 m^4 M^6-11 m^2 M^8+4 M^10))/(3 M^10 (-2 m^2+M^2)^3));
 
 g21[m_,M_]:=-((16 (-8 m^8+12 m^6 M^2-8 m^4 M^4+3 m^2 M^6+M^8))/(3 M^10 (-2 m^2+M^2)^2));
 
 
 
 (g32[m,M]-g33[m,M]*3)/(2*g22[m,M]-g21[m,M])//FullSimplify
 
 
 (* expected result *)
 ((2J(J+1)*(M^2-2m^2)/(M^2-4m^2)-3)/(M^2-2m^2)^4)/(2/(M^2-2m^2)^3)/.{J->2}//FullSimplify


(* exchange spinless particle *)
Ms0[s_, u_, M_, m_] := g^2*M^2*( 1/(t-M^2)+1/(s-M^2)+1/(u-M^2) )/.{t -> 4m^2-s-u};

Ms0[s,u,M,m]//FullSimplify


Print["g[3,3]= ", SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m],{s, 2m^2, 0}], {u, 2m^2, 3}]//FullSimplify];

Print["g[3,2]= ", SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m],{s,2m^2,1}],{u,2m^2,2}]//FullSimplify ];

Print["g[2,1]= ", SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m],{s,2m^2,1}],{u,2m^2,1}]//FullSimplify ];

Print["g[2,2]= ", SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m],{s,2m^2,0}],{u,2m^2,2}]//FullSimplify ];


g3[M_,m_] := -3*g^2 M^2 (1/M^8-1/(-2 m^2+M^2)^4)+ (3 g^2)/M^6;

g2[M_,m_] := (2 g^2)/M^4+2*g^2 M^2 (-(1/M^6)+1/(2 m^2-M^2)^3)//FullSimplify

g2[M,m]//FullSimplify

g3[M,m]/g2[M,m]//FullSimplify

Plot[g3[1,m]/g2[1,m],{m,0,0.45},PlotTheme->"Detailed"]





