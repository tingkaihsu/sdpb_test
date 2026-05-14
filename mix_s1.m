(* ::Package:: *)

Ms2[s_, u_, M_, m_] := ((s-u)^2-1/3 (4m^2-t)^2)/(t-M^2)+((t-u)^2-1/3 (4m^2-s)^2)/(s-M^2)+((s-t)^2-1/3 (4m^2-u)^2)/(u-M^2)/.{t -> 4m^2-s-u};


Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m]/.{M -> M/e}, {e, 0, 6}], {s, 2m^2, 1}], {u, 2m^2, 2}]//FullSimplify ];

Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[SeriesCoefficient[Ms2[s,u,M,m]/.{M->M/e},{e,0,4}],{s,2m^2,0}],{u,2m^2,2}]//FullSimplify ];


 -(-((32 m^2)/M^6)/(-((52 m^2)/(3 M^4))))*M^2//N


Ms0[s_, u_, M_, m_] := 1/(t-M^2)+1/(s-M^2)+1/(u-M^2)/.{t -> 4m^2-s-u};

Ms0[s,u,M,m]//FullSimplify


Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m]/.{M -> M/e}, {e, 0, 8}],{s,2m^2,2}],{u,2m^2,1}]//FullSimplify ];

Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[SeriesCoefficient[Ms0[s,u,M,m]/.{M->M/e},{e, 0, 6}],{s,2m^2,1}],{u,2m^2,1}]//FullSimplify ];


-(3/M^8)/( -(2/M^6))*M^2//N//FullSimplify
