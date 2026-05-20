maVal = SetPrecision[0.150, 1700];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 1700];
  mA = N[maVal, 1700];
  N[(81*Sqrt[sp/(-4*mA^2 + sp)]*(8*mA^4*(14*mA^2 - 15*sp)*(-8*mA^2 + 3*sp)^(3/2)*Hypergeometric2F1[-J, 1 + J, 1, (4*mA^2)/(12*mA^2 - 3*sp)] + (8*mA^4 - 18*mA^2*sp + 9*sp^2)*((-2*I)*mA*(10*mA^2 - 9*sp)*(8*mA^2 - 3*sp)*LegendreP[J, 1, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))] + Sqrt[-8*mA^2 + 3*sp]*(-8*mA^4 + 18*mA^2*sp - 9*sp^2)*LegendreP[J, 2, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))])))/(4*mA^2*(-2*mA^2 + sp)*(-8*mA^2 + 3*sp)^(3/2)*(8*mA^4 - 18*mA^2*sp + 9*sp^2)^3), 1700]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 1700];
  mA = N[maVal, 1700];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 1700]
];

X53[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 1700];
  mA = N[maVal, 1700];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*((-4 + J)*(-2 + J)*(3 + J)*(5 + J) + ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 1700]
];

Print["n4[J = 10000] = ", n4[0.000001`1700, 10000]//FullSimplify]

Print["X52[J = 10000] = ", X52[0.000001`1700, 10000]//FullSimplify]

Print["X53[J = 10000] = ", X53[0.000001`1700, 10000]//FullSimplify]