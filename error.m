G4 = g4*(s^2+t^2+u^2)^2/.{u -> -s-t};

Print["t^4 coefficient = ", SeriesCoefficient[G4, {t,0,4}]//FullSimplify];

Print["s*t^3 coefficient = ", SeriesCoefficient[SeriesCoefficient[G4, {t,0,3}], {s,0,1}]//FullSimplify];

Print["s^2*t^2 coefficient = ", SeriesCoefficient[SeriesCoefficient[G4, {t,0,2}], {s,0,2}]//FullSimplify];

Print[J*(J+1)*(2J*(J+1)-(5d-4))/.{d->4}//FullSimplify]

G3 = g3(s*t*u)/.{u -> -s-t};

Print["t^3 coefficient = ", SeriesCoefficient[G3, {t,0,3}]//FullSimplify];

Print["s*t^2 coefficient = ", SeriesCoefficient[SeriesCoefficient[G3, {t,0,2}], {s,0,1}]//FullSimplify];

(* forward dispersive sum rule *)
KerfwdMsum[sp_, t_, mA_, s1_, s2_, k_, J_] := ( 1/((sp-2mA^2)*((sp-s1)*(sp-s2))^(k/2)) - 1/(((4mA^2-sp-t)-2mA^2)*(((4mA^2-sp-t)-s1)*((4mA^2-sp-t)-s2))^(k/2)) ) * Sqrt[sp/(sp-4mA^2)] * LegendreP[J, 1+2*t/(sp-4mA^2)];

Print[""]

Print["at fwd expansion points, g[2,0] = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, 2mA^2, 2mA^2-t, 2, 8], {sp, Infinity}], {t, 0, 0}]//FullSimplify]
Print["(twice-subtracted) sum rule at fwd = ", SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 0}]//FullSimplify]

G2p = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 0}]//FullSimplify;

Print["Above under the massless limit = ", Limit[G2p, mA->0]//FullSimplify]

Print[""]

Print["at fwd expansion points, g[3,1] = ", SeriesCoefficient[-Residue[KerfwdMlow[sp, s, t, mA, 2mA^2, 2mA^2-t, 2, 8], {sp, Infinity}], {t, 0, 1}]//FullSimplify]
G3p = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 1}]//FullSimplify;
Print["(twice-subtracted) sum rule at fwd = ", SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 1}]//FullSimplify]

Print["Above under the massless limit = ", Limit[G3p, mA->0]//FullSimplify]

Print[""]

n4 = SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 2, J], {t, 0, 2}] - 2*SeriesCoefficient[KerfwdMsum[sp, t, mA, 2mA^2, 2mA^2-t, 4, J], {t, 0, 0}]//FullSimplify;

Print["n4 = ", n4]

Print["large J limit of n4 = ", SeriesCoefficient[n4, {J, 0, 4}] // Normal//FullSimplify];

Print[""]