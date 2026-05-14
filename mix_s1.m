(* ::Package:: *)

Ms2[s_, u_, M_, m_, g_:1] := Module[{t = 4 m^2 - s - u},
  g^2 (
    (3 (s - u)^2 - (4 m^2 - t)^2)/(12 (t - M^2)) +
    (3 (t - u)^2 - (4 m^2 - s)^2)/(12 (s - M^2)) +
    (3 (s - t)^2 - (4 m^2 - u)^2)/(12 (u - M^2))
  )
]


Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1], {s, 2m^2, 1}], {u, 2m^2, 2}]//FullSimplify ];

Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m,1],{s, 2m^2, 1}],{u, 2m^2, 1}]//FullSimplify ];


 g3[m_,M_] := -((4 m^4)/M^8)+(2 m^2)/M^6-1/(2 M^4)+M^2/(2 m^2-M^2)^3-1/(-2 m^2+M^2)^2;
 
 g2[m_,M_]:= (2 (4 m^4-2 m^2 M^2+M^4-(3 M^8)/(-2 m^2+M^2)^2))/(3 M^6);
 
 Plot[-(g3[m,1]/g2[m,1]),{m,0,0.45}, PlotTheme->"Detailed"]
 
 Solve[D[-(g3[m,1]/g2[m,1]),m]==0,{m},Assumptions->{m>0}]//N//FullSimplify


Ms0[s_, u_, M_, m_] := 1/(t-M^2)+1/(s-M^2)+1/(u-M^2)/.{t -> 4m^2-s-u};

Ms0[s,u,M,m]//FullSimplify


Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m],{s,2m^2,1}],{u,2m^2,2}]//FullSimplify ];

Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m],{s,2m^2,1}],{u,2m^2,1}]//FullSimplify ];


-(3/M^8)/( -(2/M^6))*M^2//N//FullSimplify
