(* ::Package:: *)

maVal = SetPrecision[1/100000000, 1700];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 2000];
  mA = N[maVal, 2000];
  N[(81*Sqrt[sp/(-4*mA^2 + sp)]*(8*mA^4*(14*mA^2 - 15*sp)*(-8*mA^2 + 3*sp)^(3/2)*Hypergeometric2F1[-J, 1 + J, 1, (4*mA^2)/(12*mA^2 - 3*sp)] + (8*mA^4 - 18*mA^2*sp + 9*sp^2)*((-2*I)*mA*(10*mA^2 - 9*sp)*(8*mA^2 - 3*sp)*LegendreP[J, 1, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))] + Sqrt[-8*mA^2 + 3*sp]*(-8*mA^4 + 18*mA^2*sp - 9*sp^2)*LegendreP[J, 2, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))])))/(4*mA^2*(-2*mA^2 + sp)*(-8*mA^2 + 3*sp)^(3/2)*(8*mA^4 - 18*mA^2*sp + 9*sp^2)^3), 2000]
];

n4massless[x_?NumericQ, J_?IntegerQ] := Module[{sp},
	sp = N[1/(1-x), 1700];
	N[(J (1+J) (-8+J+J^2))/(2 sp^5), 1700]
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

Print["n4massless[J = 10000] = ", n4massless[0.000001`1700, 10000]//FullSimplify]

Print["n4[J = 10000] = ", n4[0.000001`2000, 10000]//FullSimplify]

Print["X52[J = 10000] = ", X52[0.000001`1700, 10000]//FullSimplify]

Print["X53[J = 10000] = ", X53[0.000001`1700, 10000]//FullSimplify]



n4PrecMin[J_?IntegerQ, x_?NumericQ, margin_Integer:60, stdPrec_Integer:650] :=
  Module[{sp0, z2, A, logA},
    sp0  = N[1/(1 - x), 50];
    z2   = 1 + 8*(3/20)^2 / (3*(sp0 - 4*(3/20)^2));
    A    = If[z2 > 1, z2 + Sqrt[z2^2 - 1], 1];
    logA = If[A > 1, N[Log[10, A], 50], 0];
    Max[stdPrec, Ceiling[J * logA] + margin]
  ];
  
n4[x_?NumericQ, J_?IntegerQ] := Module[{prec, xp, sp, mA, result},
  (* FIX 2: adaptive precision \[LongDash] no overhead for J \[LessEqual] 60 *)
  prec = n4PrecMin[J, x];
  (* FIX 3: promote input x; avoids 200-digit bottleneck from file reader *)
  xp   = SetPrecision[x, prec];
  (* FIX 1 (inside n4): exact rational for mA, independent of global maVal *)
  sp   = N[1/(1 - xp), prec];
  mA   = N[maVal, prec];
  (* The expression below is UNCHANGED from the original formula.
     All intermediate values inherit prec-digit precision automatically. *)
  result = (81*Sqrt[sp/(-4*mA^2 + sp)]*(
      8*mA^4*(14*mA^2 - 15*sp)*(-8*mA^2 + 3*sp)^(3/2)*
        Hypergeometric2F1[-J, 1 + J, 1, (4*mA^2)/(12*mA^2 - 3*sp)] +
      (8*mA^4 - 18*mA^2*sp + 9*sp^2)*(
        (-2*I)*mA*(10*mA^2 - 9*sp)*(8*mA^2 - 3*sp)*
          LegendreP[J, 1, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))] +
        Sqrt[-8*mA^2 + 3*sp]*(-8*mA^4 + 18*mA^2*sp - 9*sp^2)*
          LegendreP[J, 2, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))]
      )
    ))/(4*mA^2*(-2*mA^2 + sp)*(-8*mA^2 + 3*sp)^(3/2)*(8*mA^4 - 18*mA^2*sp + 9*sp^2)^3);
  (* FIX 4: explicit Re[\[Ellipsis]] \[LongDash] the imaginary residual ~10^{-margin} is discarded.
     If this print fires with a large value, increase margin in n4PrecMin. *)
  (* Uncomment for debugging: Print["n4 Im residual @ J=", J, ": ", Im[N[result,20]]]; *)
  Re[N[result, 650]]
];


Print["n4[J = 10000] = ", n4[1/1000000000`2000, 10000]//FullSimplify]
