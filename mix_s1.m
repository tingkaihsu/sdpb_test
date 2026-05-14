Mfull[s_, u_, M_, m_] := (u-t)/(s-M^2) + (s-u)/(t-M^2) + (t-s)/(u-M^2)/.{t -> 4m^2-s-u};


Print["Large-M expansion: ", SeriesCoefficient[SeriesCoefficient[SeriesCoefficient[Mfull[s,u,M,m]/.{M -> M/e}, {e, 0, 6}], {s, 2m^2, 1}], {u, 2m^2, 2}]//FullSimplify ];

Print["Large-M expansion: ", SeriesCoefficient[Mfull[s,u,M,m]/.{M -> M/e}, {e, 0, 2}]//FullSimplify ];