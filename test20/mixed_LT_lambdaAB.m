(* ::Package:: *)

(* ======================================================================== *)
(* mixed_LT_lambdaAB.m                                                      *)
(*                                                                          *)
(* Linear-trajectory bootstrap of the mixed A (massive, Z2 +) / B           *)
(* (massless, Z2 -) system, D = 10.  Bound on                               *)
(*                                                                          *)
(*      lambda_AA lambda_BB / Sqrt[ghatA ghatB]   at fixed r = ghatA/ghatB  *)
(*                                                                          *)
(* for the isolated X2 (m2sq, J2).  Notation: mixed_state_bootstrap_v2.tex. *)
(*                                                                          *)
(*   v   = {ghatB, ghatA, lambda22, nAAAA, nBBBB, nAB}                      *)
(*   a   = (a1, a2, -1, alpha, beta, gamma)          (upper bound)          *)
(*   a.v >= 0 (PSD) on every state,  maximize a.obj = -(a1 + a2 r)          *)
(*   SDPB objective V  =>  lambda lambda/Sqrt[ghatA ghatB] <= -V/Sqrt[r]    *)
(*                                                                          *)
(* Usage:  wolframscript -file mixed_LT_lambdaAB.m 1/5 1                    *)
(*         (mAsq, r; exact numbers.  Defaults 1/5 and 1.)                   *)
(* ======================================================================== *)

LTDirectory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{LTDirectory, "..", "SDPB.m"}]];


(* ------------------------------------------------------------------------ *)
(* BLOCK 1.  Input                                                          *)
(* ------------------------------------------------------------------------ *)

mAsq = If[Length[$ScriptCommandLine] >= 2,
  ToExpression[StringTrim[$ScriptCommandLine[[2]]]], 1/5];
rRatio = If[Length[$ScriptCommandLine] >= 3,
  ToExpression[StringTrim[$ScriptCommandLine[[3]]]], 1];

mA     = Sqrt[mAsq];
JA     = 0;
m2sq   = 1;        (* X2, units m2 = 1                                    *)
J2     = 2;
mgapSq = 2;        (* neutral continuum M^2 >= mgapSq, J even             *)
mcSq   = 2;        (* charged continuum M^2 >= mcSq,   all J              *)

boundType = "upper";   (* "upper" or "lower"                              *)
includeAB = True;      (* False = gamma = 0 regression test (Sec. 5.5)    *)

If[!(MatchQ[mAsq, _Integer | _Rational] && MatchQ[rRatio, _Integer | _Rational] &&
     0 < mAsq && 4 mAsq < m2sq < mgapSq && mAsq < mcSq && rRatio > 0),
  Print["Invalid input: need exact 0 < 4 mAsq < m2sq < mgapSq, mAsq < mcSq, r > 0."];
  Exit[2]
];
Print["mAsq = ", mAsq, ", r = ", rRatio, ", bound = ", boundType,
      ", includeAB = ", includeAB];


(* ------------------------------------------------------------------------ *)
(* BLOCK 2.  Null indices                                                   *)
(*                                                                          *)
(* row n  <->  (k,q), the coefficient of (s-s0)^(k-q) t^q                   *)
(*   AAAA  (s0 = 0)      n = 0..3 : (1,0) (2,0) (3,0) (3,1)  = test20       *)
(*   BBBB  (s0 = 0)      n = 0..2 : (1,0) (2,0) (3,0)                       *)
(*         [(3,1) vanishes identically for massless externals]              *)
(*   AB    (s0 = -mA^2)  n = 0..5 : (0,0) (1,0) (1,1) (2,0) (2,1) (2,2)     *)
(* ------------------------------------------------------------------------ *)

aaNullIndices = Range[0, 3];
bbNullIndices = Range[0, 2];
abNullIndices = If[includeAB, Range[0, 5], {}];

nAA = Length[aaNullIndices];
nBB = Length[bbNullIndices];
nAB = Length[abNullIndices];

zeroAA = ConstantArray[0, nAA];
zeroBB = ConstantArray[0, nBB];
zeroAB = ConstantArray[0, nAB];

functionalDimension = 3 + nAA + nBB + nAB;  (* {ghatB, ghatA, lambda22, nulls} *)


(* ------------------------------------------------------------------------ *)
(* BLOCK 3.  Gegenbauer helpers, D = 10 (lambda = 7/2), P_J(1) = 1          *)
(* ------------------------------------------------------------------------ *)

(* p_n = P_J^(n)(1): polynomial in J *)
dP1[n_, J_] := Product[(J - i) (J + 7 + i), {i, 0, n - 1}]/(2^n Pochhammer[4, n]);

(* P_J^(n)(z) at a general point, J an integer *)
dPz[n_, J_, z_] := If[n > J, 0,
  2^n Pochhammer[7/2, n] GegenbauerC[J - n, 7/2 + n, z]/GegenbauerC[J, 7/2, 1]];

(* AB kinematics at the null point (s,t) = (-mA^2, 0) *)
alphaAB[x_] := 2 x/(x - mAsq)^2;              (* z_AB = 1 + alphaAB t          *)
zc[x_]      := 1 - alphaAB[x] (x - 3 mAsq);   (* charged u'-channel angle      *)
zn[x_]      := Sqrt[(x - 4 mAsq)/x];          (* neutral AA -> BB angle        *)
rN[x_]      := Sqrt[x (x - 4 mAsq)];


(* ------------------------------------------------------------------------ *)
(* BLOCK 4.  AA -> AA rows   (KX, NA of AAAA_test20.m, verbatim)            *)
(* ------------------------------------------------------------------------ *)

KAAAA[n_, x_, J_] := {(32*mA^4 - 16*mA^2*x + 14*J*mA^2*x + 2*J^2*mA^2*x + 2*x^2 - 7*J*x^2 - J^2*x^2)/(2*x^2*(-4*mA^2 + x)^2), (-1280*mA^6 + 960*mA^4*x - 240*mA^2*x^2 - 112*J*mA^2*x^2 + 82*J^2*mA^2*x^2 + 28*J^3*mA^2*x^2 + 2*J^4*mA^2*x^2 + 20*x^3 + 91*J*x^3 - 36*J^2*x^3 - 14*J^3*x^3 - J^4*x^3)/(20*x^3*(-4*mA^2 + x)^3), (92160*mA^8 - 92160*mA^6*x + 34560*mA^4*x^2 - 5760*mA^2*x^3 + 2016*J*mA^2*x^3 - 2260*J^2*mA^2*x^3 - 42*J^3*mA^2*x^3 + 242*J^4*mA^2*x^3 + 42*J^5*mA^2*x^3 + 2*J^6*mA^2*x^3 + 360*x^4 - 2142*J*x^4 + 1409*J^2*x^4 + 147*J^3*x^4 - 112*J^4*x^4 - 21*J^5*x^4 - J^6*x^4)/(360*x^4*(-4*mA^2 + x)^4), -1/5*(J*(7 + J)*mA^2*(80*mA^4 - 76*mA^2*x + 14*J*mA^2*x + 2*J^2*mA^2*x + 23*x^2 - 7*J*x^2 - J^2*x^2))/((4*mA^2 - x)^4*x^3)}[[n+1]];

NAAAA[n_, x_, J_] := {1/mA^4, 1/mA^6, 1/mA^8, 0}[[n+1]];


(* ------------------------------------------------------------------------ *)
(* BLOCK 5.  BB -> BB rows                                                  *)
(* ------------------------------------------------------------------------ *)

KBBBB[n_, x_, J_] := {
  (1 - 4 dP1[1, J])/x^2,
  (1 + 2 dP1[1, J] - 4 dP1[2, J])/x^3,
  (1 - 2 dP1[1, J] + 2 dP1[2, J] - 8/3 dP1[3, J])/x^4
}[[n+1]];

NBBBB[n_, x_, J_] := {1/mA^4, 1/mA^6, 1/mA^8}[[n+1]];   (* A exchange, gABB^2 *)


(* ------------------------------------------------------------------------ *)
(* BLOCK 6.  AB -> AB rows, charged states (rho_AB^2)                       *)
(*   X = M^2 + mA^2,  w = M^2 - 3 mA^2,  p_n = P^(n)(1),  c_n = P^(n)(zc)   *)
(*   The 1/w poles are removable (they cancel between R1 and R2).           *)
(* ------------------------------------------------------------------------ *)

KABAB[n_, x_, J_] := With[
  {X = x + mAsq, w = x - 3 mAsq, al = alphaAB[x],
   p1 = dP1[1, J], p2 = dP1[2, J],
   c0 = dPz[0, J, zc[x]], c1 = dPz[1, J, zc[x]], c2 = dPz[2, J, zc[x]]},
  {
    1/X + (1 - c0)/w,
    1/X^2 - (1 - c0)/w^2 + al c1/w,
    al p1 (1/X + 1/w) - (1 - c0)/w^2,
    1/X^3 + (1 - c0)/w^3 - al c1/w^2 - al^2 c2/(2 w),
    2 (1 - c0)/w^3 + al p1 (1/X^2 - 1/w^2) - al c1/w^2,
    (1 - c0)/w^3 - al p1/w^2 + al^2 p2/2 (1/X + 1/w)
  }[[n+1]]
];

NABAB[n_, x_, J_] := {1/mA^2, 1/mA^4, 0, 1/mA^6, 0, 0}[[n+1]];   (* B exchange, gABB^2 *)


(* ------------------------------------------------------------------------ *)
(* BLOCK 7.  AB -> AB rows, neutral states (rho_AA rho_BB)                  *)
(*   e_n = P^(n)(zn),  rN = Sqrt[M^2 (M^2 - 4 mA^2)]                        *)
(* ------------------------------------------------------------------------ *)

KAABB[n_, x_, J_] := With[
  {r = rN[x], e0 = dPz[0, J, zn[x]], e1 = dPz[1, J, zn[x]], e2 = dPz[2, J, zn[x]]},
  {
    -e0/x,
    -2 e1/(r x),
    -e0/x^2,
    -2 e2/(r^2 x),
    -2 e1/(r x^2),
    -e0/x^3
  }[[n+1]]
];

NAABB[n_, x_, J_] := {-1/mA^2, 0, -1/mA^4, 0, 0, -1/mA^6}[[n+1]];   (* A, gAAA gABB *)


(* ------------------------------------------------------------------------ *)
(* BLOCK 8.  Light rows = heavy rows at the light mass (consistency check)  *)
(*   A = (mAsq, 0) in AAAA, BBBB, AABB;  B = (0, 0) in ABAB.                *)
(* ------------------------------------------------------------------------ *)

lightRowsOK[] := And @@ Flatten[{
  Table[PossibleZeroQ[KAAAA[n, mAsq, 0] - NAAAA[n, mAsq, 0]], {n, aaNullIndices}],
  Table[PossibleZeroQ[KBBBB[n, mAsq, 0] - NBBBB[n, mAsq, 0]], {n, bbNullIndices}],
  Table[PossibleZeroQ[KABAB[n, 0, 0]    - NABAB[n, 0, 0]],    {n, abNullIndices}],
  Table[PossibleZeroQ[KAABB[n, mAsq, 0] - NAABB[n, mAsq, 0]], {n, abNullIndices}]
}];

If[!lightRowsOK[], Print["Light rows do not match heavy rows at the light mass."]; Exit[3]];
Print["light-row check passed"];


(* ------------------------------------------------------------------------ *)
(* BLOCK 9.  Normalization and phase-space weights                          *)
(* ------------------------------------------------------------------------ *)

g0AAWeight[x_] := 2/(x - 2 mAsq);     (* ghatA at (s,t) = (2mA^2, 0) *)
g0BBWeight[x_] := 2/x;                (* ghatB at (s,t) = (0, 0)     *)

(* phi = (2p)^7/Sqrt[x].  Neutral element ij carries Sqrt[phi_i phi_j]:   *)
(* congruence diag(Sqrt[phiAA], Sqrt[phiBB]), PSD cone unchanged.         *)
phiAA[x_] := (x - 4 mAsq)^(7/2)/Sqrt[x];
phiBB[x_] := x^3;
phiAB[x_] := (x - mAsq)^7/x^4;

spinScale[j_] := 1;


(* ------------------------------------------------------------------------ *)
(* BLOCK 10.  Component vectors {ghatB, ghatA, lambda22, nAAAA, nBBBB, nAB} *)
(* ------------------------------------------------------------------------ *)

aaRow[x_, j_]   := Table[KAAAA[n, x, j], {n, aaNullIndices}];
bbRow[x_, j_]   := Table[KBBBB[n, x, j], {n, bbNullIndices}];
abRow[x_, j_]   := Table[KABAB[n, x, j], {n, abNullIndices}];
aabbRow[x_, j_] := Table[KAABB[n, x, j], {n, abNullIndices}];

(* neutral state; lam12 = entry of E = {{0,1/2},{1/2,0}}: 1/2 for X2, else 0 *)
vecAA[x_, j_]          := Join[{0, g0AAWeight[x], 0}, aaRow[x, j], zeroBB, zeroAB];
vecBB[x_, j_]          := Join[{g0BBWeight[x], 0, 0}, zeroAA, bbRow[x, j], zeroAB];
vecAABB[x_, j_, lam12_] := Join[{0, 0, lam12}, zeroAA, zeroBB, aabbRow[x, j]/2];

(* charged state *)
vecABAB[x_, j_] := Join[{0, 0, 0}, zeroAA, zeroBB, abRow[x, j]];

(* light A and B, acting on (gAAA, gABB); A is absorbed into the hats *)
vecLightAA = Join[{0, 0, 0}, Table[NAAAA[n, mAsq, JA], {n, aaNullIndices}], zeroBB, zeroAB];
vecLightBB = Join[{0, 0, 0}, zeroAA, Table[NBBBB[n, mAsq, JA], {n, bbNullIndices}],
                  Table[NABAB[n, 0, 0], {n, abNullIndices}]];
vecLightAB = Join[{0, 0, 0}, zeroAA, zeroBB, Table[NAABB[n, mAsq, JA]/2, {n, abNullIndices}]];


(* ------------------------------------------------------------------------ *)
(* BLOCK 11.  PMP blocks = the four inequalities of Sec. 5.4                *)
(*   PolyA  (i)   light A and B         2x2 constant                        *)
(*   PolyH  (ii)  X2                    2x2,  F - E >= 0                    *)
(*   PolyX  (iii) neutral continuum     2x2,  J even                        *)
(*   PolyC  (iv)  charged continuum     1x1,  all J                         *)
(* ------------------------------------------------------------------------ *)

prefactor[y_] := DampedRational[1, {}, 1/E, y];

neutralMatrix[j_, x_, lam12_] := Module[{wAA, wBB, off},
  wAA = phiAA[x];
  wBB = phiBB[x];
  off = Sqrt[wAA wBB] vecAABB[x, j, lam12];
  spinScale[j] {{wAA vecAA[x, j], off}, {off, wBB vecBB[x, j]}}
];

PolyA[y_] := PositiveMatrixWithPrefactor[
  prefactor[y],
  {{vecLightAA, vecLightAB}, {vecLightAB, vecLightBB}}
];

PolyH[j_, x_, y_] := PositiveMatrixWithPrefactor[prefactor[y], neutralMatrix[j, x, 1/2]];

PolyX[j_, x_, y_] := PositiveMatrixWithPrefactor[prefactor[y], neutralMatrix[j, x, 0]];

PolyC[j_, x_, y_] := PositiveMatrixWithPrefactor[
  prefactor[y],
  {{spinScale[j] phiAB[x] vecABAB[x, j]}}
];


(* ------------------------------------------------------------------------ *)
(* BLOCK 12.  Assemble and write                                            *)
(* ------------------------------------------------------------------------ *)

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
  {
    npts = 150, Jmax = 100, sgn,
    phiSamples, massSamples, evenSpins, allSpins, chargedX,
    pA, pH, pX, pC, pols, norm, obj
  },

  phiSamples = N[Table[Pi/2 + Pi/2 Cos[(k + 1/2) Pi/npts], {k, 0, npts - 1}], prec];
  massSamples = Sin[#/2]^2 & /@ phiSamples;
  evenSpins = Range[0, Jmax, 2];
  allSpins  = Range[0, Jmax];

  (* charged samples: stay off the removable point M^2 = 3 mA^2 *)
  chargedX = Select[mcSq/(1 - #) & /@ massSamples, Abs[# - 3 mAsq] > 10^-10 &];

  pA = {N[PolyA[x], prec]};
  pH = {N[PolyH[J2, m2sq, x], prec]};
  pX = Flatten[N[ParallelTable[PolyX[i, mgapSq/(1 - m), x],
         {i, evenSpins}, {m, massSamples}], prec], 1];
  pC = If[nAB == 0, {},
         Flatten[N[ParallelTable[PolyC[i, xc, x],
           {i, allSpins}, {xc, chargedX}], prec], 1]];

  pols = Join[pA, pH, pX, pC];

  (* upper: a_3 = -1 ;  lower: a_3 = +1 (E -> -E) *)
  sgn  = If[boundType === "upper", 1, -1];
  norm = -sgn N[Join[{0, 0, 1}, zeroAA, zeroBB, zeroAB], prec];
  obj  = -N[Join[{1, rRatio, 0}, zeroAA, zeroBB, zeroAB], prec];

  If[
    Length[norm] =!= functionalDimension ||
    Length[obj]  =!= functionalDimension,
    Print["Objective or normalization dimension mismatch."];
    Abort[]
  ];

  Print["functional dimension = ", functionalDimension,
        "  (3 + ", nAA, " AAAA + ", nBB, " BBBB + ", nAB, " AB)"];
  Print["number of PMP blocks = ", Length[pols],
        "  (light 1, X2 1, neutral ", Length[pX], ", charged ", Length[pC], ")"];
  Print["bound = ", If[sgn == 1, "-V/Sqrt[r]", "+V/Sqrt[r]"],
        " with V the SDPB objective; Sqrt[r] = ", N[Sqrt[rRatio], 10]];

  WritePmpJson[
    datfile,
    SDP[obj, norm, pols],
    prec,
    getAnalyticSampleData
  ]
];

PMP2SDP[FileNameJoin[{LTDirectory,"n_pmp.json"}], 1000];
