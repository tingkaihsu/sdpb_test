f[sp_, s_, t_, mA_] := ( 1/(sp - s) + 1/(sp - u) )*LegendreP[J, 1+2*t/(sp-4mA^2)]/.{u -> 4mA^2 - s - t}

ma = 0;

g20[x_, J_] := 1/2*Sqrt[ sp/ (sp-4*mA^2) ] * (sp^(-3) + (-4*mA^2 + sp)^(-3))/.{sp -> 1/(1-x), mA -> ma};

g31[x_, J_] := -Sqrt[ sp/ (sp-4*mA^2) ] * ((-3 + J*(1 + J)*(-4*mA^2 + sp)^3*(sp^(-3) + (-4*mA^2 + sp)^(-3)))/(-4*mA^2 + sp)^4)/.{sp -> 1/(1-x), mA -> ma};

Print["g[2, 0] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, mA], {s, 0, 2}], {t, 0, 0}]//FullSimplify];

Print["massless limit g[2, 0] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, 0], {s, 0, 2}], {t, 0, 0}]//FullSimplify];

Print["g[3, 1] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, mA], {s, 0, 2}], {t, 0, 1}]//FullSimplify];

Print["massless limit g[3, 1] = ", SeriesCoefficient[SeriesCoefficient[f[sp, s, t, 0], {s, 0, 2}], {t, 0, 1}]//FullSimplify];

Print["massless limit g[2, 0] = ", g20[x, J]//FullSimplify]

Print["massless limit g[3, 1] = ", g31[x, J]//FullSimplify]