int = Integrate[ Exp[I*k*r*x], {x, -1, 1}];


Print[ FullSimplify[ int ] ]

int = Integrate[ Sin[k*r] / k * V0 * r, {r, 0, a}, Assumptions -> {k > 0, V0 > 0, a > 0}] // FullSimplify;

Print[ int ]

f = (V0*(-(a*k*Cos[a*k]) + Sin[a*k]))/k^3 * 2*m / h^2 /. {k -> 2*k*Sin[\[Theta]/2]};


Print[ FullSimplify[ SeriesCoefficient[f, {k, 0, 0}] ] ] 
