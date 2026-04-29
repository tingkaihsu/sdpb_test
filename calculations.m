f[sp_, s_, t_, mA_] := ( 1/(sp - s) + 1/(sp - u) )*LegendreP[J, 1+2*t/(sp-4mA^2)]/.{u -> 4mA^2 - s - t}

Print["g[2, 0] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, mA], {s, 0, 2}], {t, 0, 0}]//FullSimplify];

Print["massless limit g[2, 0] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, 0], {s, 0, 2}], {t, 0, 0}]//FullSimplify];

Print["g[3, 1] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, mA], {s, 0, 2}], {t, 0, 1}]//FullSimplify];

Print["massless limit g[3, 1] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, 0], {s, 0, 2}], {t, 0, 1}]//FullSimplify];