(* ::Package:: *)

(* pmp generator for mixed scattering *)
(* bootstrapping the Wilson-coefficient island with single-channel null constraints *)

ClearAll[NumericalPositiveMatrixWithPrefactor];

toJsonNestedNumberArray[expr_, prec_] := expr /. n_?NumericQ :> toJsonNumber[n, prec];

toJsonObject[NumericalPositiveMatrixWithPrefactor[pmp_?AssociationQ], prec_, getSampleDataFn_:Function[<||>]] :=
Module[
  {sampleData, samplePoints, sampleScalings, functionValues, basisValues,
   prefactor, reducedPrefactor},

  sampleData = getSampleDataFn[NumericalPositiveMatrixWithPrefactor[pmp], prec];

  prefactor        = Lookup[pmp, "prefactor",        1];
  reducedPrefactor = Lookup[pmp, "reducedPrefactor", prefactor];

  samplePoints   = Lookup[sampleData, "samplePoints",  Lookup[pmp, "samplePoints",  Missing[]]];
  sampleScalings = Lookup[sampleData, "sampleScalings", Lookup[pmp, "sampleScalings", Missing[]]];
  functionValues = Lookup[pmp, "polynomials", Missing[]];
  basisValues    = Lookup[sampleData, "basisValues", Lookup[pmp, "basisValues", Missing[]]];

  DeleteMissing @ <|
    "prefactor"        -> toJsonDampedRational[prefactor,        prec],
    "reducedPrefactor" -> toJsonDampedRational[reducedPrefactor, prec],
    "samplePoints"     -> toJsonNumberArray[samplePoints,   prec],
    "sampleScalings"   -> toJsonNumberArray[sampleScalings, prec],
    "polynomials"      -> If[MissingQ[functionValues], Missing[],
      toJsonNestedNumberArray[functionValues, prec]
    ],
    "basisValues"      -> If[MissingQ[basisValues], Missing[], toJsonNestedNumberArray[basisValues, prec]]
  |>
];

WritePmpJsonNumerical[
  file_,
  SDP[objective_, normalization_, positiveMatricesWithPrefactors_],
  prec_,
  getSampleDataFn_:Function[<||>]
] := exportJson[
  file,
  <|
    "objective"     -> toJsonNumberArray[objective,     prec],
    "normalization" -> toJsonNumberArray[normalization, prec],
    "PositiveMatrixWithPrefactorArray" ->
      Table[toJsonObject[pmp, prec, getSampleDataFn], {pmp, positiveMatricesWithPrefactors}]
  |>
];

<< "../SDPB.m";

(* problem-specific *)
(* let the mass be 4mA^2 < M^2 = 1 where M = 1 is the first isolated massive pole, and the mass gap is Mgap = 2 *)
(* 
   LINEAR TRAJECTORY CONDITION (from arXiv:2510.07991, Eq. 13):
   For a 2-state system with (m1^2, J1) and (m2^2, J2=J1+2), the next state on
   a linear trajectory sits at m3^2 = 2*m2^2 - m1^2. For this to be in the UV
   (above M^2), we need: m1^2 <= 2*m2^2 - M^2.
   
   Current: m1^2=1, m2^2=6/5, M^2=2  =>  m1,c^2 = 2(6/5)-2 = 2/5 = 0.4
   Since m1^2=1 > 0.4, the linear trajectory CANNOT form with these parameters.
   
   To find a linear trajectory, choose one of:
     (A) Lower mgap: M^2 <= 2*m2^2 - m1^2 = 7/5  (e.g., mgap = 7/5 or 13/10)
     (B) Raise m2:   m2^2 >= (m1^2+M^2)/2 = 3/2   (e.g., m2 = 3/2 or 8/5)
     (C) Use a 3+ state ansatz including ell=6 at m3^2=7/5, etc.
*)
Print["Mass scales are normalized by the second isolated state..."]
Print[""]

J1 = 2;
m1 = N[1/2, 600];
mgap = N[2, 600];
mAval = N[1/1000, 600];


Print["J1 = ", J1];
Print["m1^2 = ", m1];
Print["m_gap^2  = ", mgap];
Print["mA = ", mAval];

(* forward limit: use our own convention *)

(* 2 delCoeff[0,2]-delCoeff[1,1] *)
g2[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((2 Sqrt[sp/(-4 mA^2+sp)])/(2 mA^2-sp)^3);
  N[result, 600]
];

(* number of null constraints = 30 *)
(* n4 *)
x42[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((243 Sqrt[-(sp/(4 mA^2-sp))] (6 (4 mA^2-sp) Sqrt[(8 mA^4-3 mA^2 sp)/(-4 mA^2+sp)^2] LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+(4 mA^2-3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (8 mA^2-3 sp)));
  N[result, 600]
];

x52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = ((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-4+J) (-2+J) (3+J) (5+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (36 mA^2+(-15+J+J^2) sp))/sp^4))/(36 (-4 mA^2+sp)^6));
  N[result, 600]
];

x62[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7);
  N[result, 600]
];

x72[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400;
  N[result, 600]
];

x73[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = ((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-6+J) (-4+J) (5+J) (7+J)-((-1+J) (2+J) (-4 mA^2+sp)^4 (64 mA^2+(-28+J+J^2) sp))/sp^5))/(576 (-4 mA^2+sp)^8);
  N[result, 600]
];

x82[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (((-8+J) (-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (9+J) (-38+J+J^2))/(4 mA^2-sp)^9+((-1+J) (2+J) (-(((-5+J) (-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) (6+J) sp^4)/(-4 mA^2+sp)^6)+129600/(-4 mA^2+sp)^2))/sp^7))/518400;
  N[result, 600]
];

x83[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = (((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-8+J) (-6+J) (-4+J) (5+J) (7+J) (9+J)+((-2+J+J^2) (-4 mA^2+sp)^4 (6400 mA^4-3200 mA^2 sp+(160-J (1+J) (-32+J+J^2)) sp^2))/sp^6))/(14400 (-4 mA^2+sp)^9));
  N[result, 600]
];

(* from x82 *)
largeJ[x_?NumericQ] := Module[{sp, mA, result},
  sp = N[mgap/(1-x), 600];
  mA = N[mAval, 600];
  result = -((Sqrt[sp/(-4 mA^2+sp)] (-32 mA^6+24 mA^4 sp-6 mA^2 sp^2+sp^3))/(259200 sp^3 (-4 mA^2+sp)^9));
  N[result, 600]
];

(* first state resonance *)
(* 2 delCoeff[0,2]-delCoeff[1,1] *)
g2fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((2 Sqrt[sp/(-4 mA^2+sp)])/(2 mA^2-sp)^3);
  N[result, 600]
];

(* number of null constraints = 30 *)
(* n4 *)
x42fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = -((243 Sqrt[-(sp/(4 mA^2-sp))] (6 (4 mA^2-sp) Sqrt[(8 mA^4-3 mA^2 sp)/(-4 mA^2+sp)^2] LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+(4 mA^2-3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (8 mA^2-3 sp)));
  N[result, 600]
];

x52fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = ((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-4+J) (-2+J) (3+J) (5+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (36 mA^2+(-15+J+J^2) sp))/sp^4))/(36 (-4 mA^2+sp)^6));
  N[result, 600]
];

x62fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7);
  N[result, 600]
];

x72fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400;
  N[result, 600]
];

x73fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = ((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-6+J) (-4+J) (5+J) (7+J)-((-1+J) (2+J) (-4 mA^2+sp)^4 (64 mA^2+(-28+J+J^2) sp))/sp^5))/(576 (-4 mA^2+sp)^8);
  N[result, 600]
];

x82fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (((-8+J) (-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (9+J) (-38+J+J^2))/(4 mA^2-sp)^9+((-1+J) (2+J) (-(((-5+J) (-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) (6+J) sp^4)/(-4 mA^2+sp)^6)+129600/(-4 mA^2+sp)^2))/sp^7))/518400;
  N[result, 600]
];

x83fst[x_?NumericQ] := Module[{sp, mA, J, result},
  sp = N[m1, 600];
  mA = N[mAval, 600];
  J = J1;
  result = (((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-8+J) (-6+J) (-4+J) (5+J) (7+J) (9+J)+((-2+J+J^2) (-4 mA^2+sp)^4 (6400 mA^4-3200 mA^2 sp+(160-J (1+J) (-32+J+J^2)) sp^2))/sp^6))/(14400 (-4 mA^2+sp)^9));
  N[result, 600]
];

M2[x_?NumericQ,J_?IntegerQ] := {
	{g2[x, J],0,0},
	{0,0,0},
	{0,0,0}
};

M2fst[x_?NumericQ] := {
  {g2fst[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M42[x_?NumericQ, J_?IntegerQ] := {
  {x42[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M42fst[x_?NumericQ] := {
  {x42fst[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M52[x_?NumericQ, J_?IntegerQ] := {
  {x52[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M52fst[x_?NumericQ] := {
  {x52fst[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M62[x_?NumericQ, J_?IntegerQ] := {
  {x62[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M62fst[x_?NumericQ] := {
  {x62fst[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M72[x_?NumericQ, J_?IntegerQ] := {
  {x72[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M72fst[x_?NumericQ] := {
  {x72fst[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M73[x_?NumericQ, J_?IntegerQ] := {
  {x73[x, J], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M73fst[x_?NumericQ] := {
  {x73fst[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

M82[x_?NumericQ, J_?IntegerQ] :={
	{x82[x, J],0,0},
	{0,0,0},
	{0,0,0}
}

M82fst[x_?NumericQ] := {
  {x82fst[x],0,0},
	{0,0,0},
	{0,0,0}
}

M83[x_?NumericQ, J_?IntegerQ] :={
	{x83[x, J],0,0},
	{0,0,0},
	{0,0,0}
}

M83fst[x_?NumericQ] := {
  {x83fst[x],0,0},
	{0,0,0},
	{0,0,0}
}

MlgeJ[x_?NumericQ] := {
  {largeJ[x], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
}

(* null *)
N0[x_?NumericQ] :={
	{0,0,0},
	{0,0,0},
	{0,0,0}
};

f11List ={
	Function[{x,J}, M2[x,J][[1,1]]],
	Function[{x,J}, M42[x,J][[1,1]]],
  Function[{x,J}, M52[x,J][[1,1]]],
  Function[{x,J}, M62[x,J][[1,1]]],
  Function[{x,J}, M72[x,J][[1,1]]],
  Function[{x,J}, M73[x, J][[1,1]]],
  Function[{x,J}, M82[x, J][[1,1]]],
  Function[{x,J}, M83[x, J][[1,1]]]
};

f22List ={
	Function[{x,J}, M2[x,J][[2,2]]],
	Function[{x,J}, M42[x,J][[2,2]]],
  Function[{x,J}, M52[x,J][[2,2]]],
  Function[{x,J}, M62[x,J][[2,2]]],
  Function[{x,J}, M72[x,J][[2,2]]],
  Function[{x,J}, M73[x,J][[2,2]]],
  Function[{x,J}, M82[x,J][[2,2]]],
  Function[{x,J}, M83[x,J][[2,2]]]
};

f33List = {
	Function[{x,J}, M2[x,J][[3,3]]],
	Function[{x,J}, M42[x,J][[3,3]]],
  Function[{x,J}, M52[x,J][[3,3]]],
  Function[{x,J}, M62[x,J][[3,3]]],
  Function[{x,J}, M72[x,J][[3,3]]],
  Function[{x,J}, M73[x,J][[3,3]]],
  Function[{x,J}, M82[x,J][[3,3]]],
  Function[{x,J}, M83[x,J][[3,3]]]
};

f12List ={
	Function[{x,J}, M2[x,J][[1,2]]],
	Function[{x,J}, M42[x,J][[1,2]]],
  Function[{x,J}, M52[x,J][[1,2]]],
  Function[{x,J}, M62[x,J][[1,2]]],
  Function[{x,J}, M72[x,J][[1,2]]],
  Function[{x,J}, M73[x,J][[1,2]]],
  Function[{x,J}, M82[x,J][[1,2]]],
  Function[{x,J}, M83[x,J][[1,2]]]
};

f21List = f12List;

f13List = {
	Function[{x,J}, M2[x,J][[1,3]]],
	Function[{x,J}, M42[x,J][[1,3]]],
  Function[{x,J}, M52[x,J][[1,3]]],
  Function[{x,J}, M62[x,J][[1,3]]],
  Function[{x,J}, M72[x,J][[1,3]]],
  Function[{x,J}, M73[x,J][[1,3]]],
  Function[{x,J}, M82[x,J][[1,3]]],
  Function[{x,J}, M83[x,J][[1,3]]]
};
f31List = f13List;

f23List = {
	Function[{x,J}, M2[x,J][[2,3]]],
	Function[{x,J}, M42[x,J][[2,3]]],
  Function[{x,J}, M52[x,J][[2,3]]],
  Function[{x,J}, M62[x,J][[2,3]]],
  Function[{x,J}, M72[x,J][[2,3]]],
  Function[{x,J}, M73[x,J][[2,3]]],
  Function[{x,J}, M82[x,J][[2,3]]],
  Function[{x,J}, M83[x,J][[2,3]]]
};

f32List = f23List;

f11FstList = {
  Function[{x}, M2fst[x][[1, 1]]],
  Function[{x}, M42fst[x][[1, 1]]],
  Function[{x}, M52fst[x][[1, 1]]],
  Function[{x}, M62fst[x][[1, 1]]],
  Function[{x}, M72fst[x][[1, 1]]],
  Function[{x}, M73fst[x][[1, 1]]],
  Function[{x}, M82fst[x][[1, 1]]],
  Function[{x}, M83fst[x][[1, 1]]]
};

f22FstList = {
  Function[{x}, M2fst[x][[2, 2]]],
  Function[{x}, M42fst[x][[2, 2]]],
  Function[{x}, M52fst[x][[2, 2]]],
  Function[{x}, M62fst[x][[2, 2]]],
  Function[{x}, M72fst[x][[2, 2]]],
  Function[{x}, M73fst[x][[2, 2]]],
  Function[{x}, M82fst[x][[2, 2]]],
  Function[{x}, M83fst[x][[2, 2]]]
};

f33FstList = {
  Function[{x}, M2fst[x][[3, 3]]],
  Function[{x}, M42fst[x][[3, 3]]],
  Function[{x}, M52fst[x][[3, 3]]],
  Function[{x}, M62fst[x][[3, 3]]],
  Function[{x}, M72fst[x][[3, 3]]],
  Function[{x}, M73fst[x][[3, 3]]],
  Function[{x}, M82fst[x][[3, 3]]],
  Function[{x}, M83fst[x][[3, 3]]]
};

f12FstList = {
  Function[{x}, M2fst[x][[1, 2]]],
  Function[{x}, M42fst[x][[1, 2]]],
  Function[{x}, M52fst[x][[1, 2]]],
  Function[{x}, M62fst[x][[1, 2]]],
  Function[{x}, M72fst[x][[1, 2]]],
  Function[{x}, M73fst[x][[1, 2]]],
  Function[{x}, M82fst[x][[1, 2]]],
  Function[{x}, M83fst[x][[1, 2]]]
};

f21FstList = f12FstList;

f13FstList = {
  Function[{x}, M2fst[x][[1, 3]]],
  Function[{x}, M42fst[x][[1, 3]]],
  Function[{x}, M52fst[x][[1, 3]]],
  Function[{x}, M62fst[x][[1, 3]]],
  Function[{x}, M72fst[x][[1, 3]]],
  Function[{x}, M73fst[x][[1, 3]]],
  Function[{x}, M82fst[x][[1, 3]]],
  Function[{x}, M83fst[x][[1, 3]]]
};

f31FstList = f13FstList;

f23FstList = {
  Function[{x}, M2fst[x][[2, 3]]],
  Function[{x}, M42fst[x][[2, 3]]],
  Function[{x}, M52fst[x][[2, 3]]],
  Function[{x}, M62fst[x][[2, 3]]],
  Function[{x}, M72fst[x][[2, 3]]],
  Function[{x}, M73fst[x][[2, 3]]],
  Function[{x}, M82fst[x][[2, 3]]],
  Function[{x}, M83fst[x][[2, 3]]]
};

f32FstList = f23FstList;

(* Large J limit *)
j11List = {
	Function[{x}, N0[x][[1, 1]]],
  Function[{x}, N0[x][[1, 1]]],
  Function[{x}, N0[x][[1, 1]]],
  Function[{x}, N0[x][[1, 1]]],
  Function[{x}, N0[x][[1, 1]]],
  Function[{x}, N0[x][[1, 1]]],
  Function[{x}, MlgeJ[x][[1, 1]]],
  Function[{x}, N0[x][[1, 1]]]
};

j22List = {
  Function[{x}, N0[x][[2, 2]]],
  Function[{x}, N0[x][[2, 2]]],
  Function[{x}, N0[x][[2, 2]]],
  Function[{x}, N0[x][[2, 2]]],
  Function[{x}, N0[x][[2, 2]]],
  Function[{x}, N0[x][[2, 2]]],
  Function[{x}, MlgeJ[x][[2, 2]]],
  Function[{x}, N0[x][[2, 2]]]
};

j33List = {
  Function[{x}, N0[x][[3, 3]]],
  Function[{x}, N0[x][[3, 3]]],
  Function[{x}, N0[x][[3, 3]]],
  Function[{x}, N0[x][[3, 3]]],
  Function[{x}, N0[x][[3, 3]]],
  Function[{x}, N0[x][[3, 3]]],
  Function[{x}, MlgeJ[x][[3, 3]]],
  Function[{x}, N0[x][[3, 3]]]
};

j12List = {
  Function[{x}, N0[x][[1, 2]]],
  Function[{x}, N0[x][[1, 2]]],
  Function[{x}, N0[x][[1, 2]]],
  Function[{x}, N0[x][[1, 2]]],
  Function[{x}, N0[x][[1, 2]]],
  Function[{x}, N0[x][[1, 2]]],
  Function[{x}, MlgeJ[x][[1, 2]]],
  Function[{x}, N0[x][[1, 2]]]
};
j21List = j12List;

j13List = {
  Function[{x}, N0[x][[1, 3]]],
  Function[{x}, N0[x][[1, 3]]],
  Function[{x}, N0[x][[1, 3]]],
  Function[{x}, N0[x][[1, 3]]],
  Function[{x}, N0[x][[1, 3]]],
  Function[{x}, N0[x][[1, 3]]],
  Function[{x}, MlgeJ[x][[1, 3]]],
  Function[{x}, N0[x][[1, 3]]]
};
j31List = j13List;

j23List = {
  Function[{x}, N0[x][[2, 3]]],
  Function[{x}, N0[x][[2, 3]]],
  Function[{x}, N0[x][[2, 3]]],
  Function[{x}, N0[x][[2, 3]]],
  Function[{x}, N0[x][[2, 3]]],
  Function[{x}, N0[x][[2, 3]]],
  Function[{x}, MlgeJ[x][[2, 3]]],
  Function[{x}, N0[x][[2, 3]]]
};
j32List = j23List;

J2 = 4;
Print["J2 = ", J2];

(* norm, which is second state resonance *)
(* 2 delCoeff[0,2]-delCoeff[1,1] *)
g2norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((2 Sqrt[sp/(-4 mA^2+sp)])/(2 mA^2-sp)^3);
  N[result, 600]
];

(* number of null constraints = 30 *)
(* n4 *)
x42norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = -((243 Sqrt[-(sp/(4 mA^2-sp))] (6 (4 mA^2-sp) Sqrt[(8 mA^4-3 mA^2 sp)/(-4 mA^2+sp)^2] LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+(4 mA^2-3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (8 mA^2-3 sp)));
  N[result, 600]
];

x52norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = ((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-4+J) (-2+J) (3+J) (5+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (36 mA^2+(-15+J+J^2) sp))/sp^4))/(36 (-4 mA^2+sp)^6));
  N[result, 600]
];

x62norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7);
  N[result, 600]
];

x72norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400;
  N[result, 600]
];

x73norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = ((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-6+J) (-4+J) (5+J) (7+J)-((-1+J) (2+J) (-4 mA^2+sp)^4 (64 mA^2+(-28+J+J^2) sp))/sp^5))/(576 (-4 mA^2+sp)^8);
  N[result, 600]
];

x82norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (J (1+J) Sqrt[sp/(-4 mA^2+sp)] (((-8+J) (-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (9+J) (-38+J+J^2))/(4 mA^2-sp)^9+((-1+J) (2+J) (-(((-5+J) (-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) (6+J) sp^4)/(-4 mA^2+sp)^6)+129600/(-4 mA^2+sp)^2))/sp^7))/518400;
  N[result, 600]
];

x83norm = Module[{sp, mA, J, result},
  sp = N[1, 600];
  mA = N[mAval, 600];
  J = J2;
  result = (((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-8+J) (-6+J) (-4+J) (5+J) (7+J) (9+J)+((-2+J+J^2) (-4 mA^2+sp)^4 (6400 mA^4-3200 mA^2 sp+(160-J (1+J) (-32+J+J^2)) sp^2))/sp^6))/(14400 (-4 mA^2+sp)^9));
  N[result, 600]
];

Jmax = 60;
Jlist = Range[0, Jmax, 2];

norm = {g2norm, x42norm, x52norm, x62norm, x72norm, x73norm, x82norm, x83norm};

obj  = {-1, 0, 0, 0, 0, 0, 0, 0};


testNumericalSDP[spFile_String, jsonFile_String, prec_:600] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, polsRegular, polsExtra},

  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " x sample points from ", spFile];
  Print["  x-points    : ", samplePoints];
  Print["  J-values    : ", Jlist, "  (", Length[Jlist], " spins, exact)"];
  Print["  large-J limit: {0,0,0,LargeJ}  (J\[Rule]\[Infinity] limit)"];

  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  polsRegular = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {
	      { 
          Table[{SetPrecision[f11List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f11List]}],
          Table[{SetPrecision[f21List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f21List]}],
          Table[{SetPrecision[f31List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f31List]}]
	      },
	      {
	        Table[{SetPrecision[f12List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f12List]}],
			    Table[{SetPrecision[f22List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f22List]}],
			    Table[{SetPrecision[f32List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f32List]}]
	      },
	      {
	        Table[{SetPrecision[f13List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f13List]}],
			    Table[{SetPrecision[f23List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f23List]}],
			    Table[{SetPrecision[f33List[[k]][samplePoints[[i]], Jlist[[j]]], prec]}, {k, Length[f33List]}]
	      }
      }
    |>],
    {i, Length[samplePoints]}, {j, Length[Jlist]}
  ];

  polsExtra = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {
	      { 
          Table[{SetPrecision[j11List[[k]][samplePoints[[i]]], prec]}, {k, Length[j11List]}],
          Table[{SetPrecision[j21List[[k]][samplePoints[[i]]], prec]}, {k, Length[j21List]}],
          Table[{SetPrecision[j31List[[k]][samplePoints[[i]]], prec]}, {k, Length[j31List]}]
	      },
	      {
	        Table[{SetPrecision[j12List[[k]][samplePoints[[i]]], prec]}, {k, Length[j12List]}],
			    Table[{SetPrecision[j22List[[k]][samplePoints[[i]]], prec]}, {k, Length[j22List]}],
			    Table[{SetPrecision[j32List[[k]][samplePoints[[i]]], prec]}, {k, Length[j32List]}]
	      },
	      {
	        Table[{SetPrecision[j13List[[k]][samplePoints[[i]]], prec]}, {k, Length[j13List]}],
			    Table[{SetPrecision[j23List[[k]][samplePoints[[i]]], prec]}, {k, Length[j23List]}],
			    Table[{SetPrecision[j33List[[k]][samplePoints[[i]]], prec]}, {k, Length[j33List]}]
	      }
      }
    |>],
    {i, Length[samplePoints]}
  ];

  polsFst = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {
	      { 
          Table[{SetPrecision[f11ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f11ShftList]}],
          Table[{SetPrecision[f21ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f21ShftList]}],
          Table[{SetPrecision[f31ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f31ShftList]}]
	      },
	      {
	        Table[{SetPrecision[f12ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f12ShftList]}],
			    Table[{SetPrecision[f22ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f22ShftList]}],
			    Table[{SetPrecision[f32ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f32ShftList]}]
	      },
	      {
	        Table[{SetPrecision[f13ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f13ShftList]}],
			    Table[{SetPrecision[f23ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f23ShftList]}],
			    Table[{SetPrecision[f33ShftList[[k]][samplePoints[[i]]], prec]}, {k, Length[f33ShftList]}]
	      }
      }
    |>],
    {i, Length[samplePoints]}
  ];

  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Join[Flatten[polsRegular], polsExtra, polsFst]],
    prec
  ];
  Print["Wrote PMP JSON to ", jsonFile]
];

Module[{myArgs, spFile, jsonFile, prec},

  myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

  If[Length[myArgs] >= 1,
    spFile   = myArgs[[1]];
    jsonFile = If[Length[myArgs] >= 2, myArgs[[2]], "numeric_pmp.json"];
    (* 600 digits exceeds SDPB 2048-bit precision (\[TildeTilde] 616.5 decimal digits)
       by a safe margin.  n4 uses higher precision adaptively when J is large. *)
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 600];

    Print["=== test12.m ==="];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    Null
  ]
];
