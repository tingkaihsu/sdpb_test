(* ::Package:: *)

(* pmp generator for mixed scattering *)

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
(* let the mass be 4mA^2 < M^2 = 1 where M  = 1 to infinity *)

maVal = N[1/10, 600];

Print["mA = ", maVal]

(* forward limit: use our own convention *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[-(2*Sqrt[sp/(-4*mA^2 + sp)])/(2*mA^2 - sp)^3, 600]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[-(Sqrt[sp/(-4*mA^2 + sp)]*(-3 - (2*J*(1 + J)*(2*mA^2 - sp))/(-4*mA^2 + sp)))/(-2*mA^2 + sp)^4, 600]
];

g4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = -((2 Sqrt[sp/(-4 mA^2+sp)])/(2 mA^2-sp)^5);
  Re[N[result, 600] ]
]

n4AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp   = N[1/(1 - x), 600];
  mA   = N[maVal, 600];
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
  Re[N[result, 600] ]
];

n4BBBB[x_?NumericQ, J_?IntegerQ] := Module[{sp, result},
  sp   = N[1/(1 - x), 600];
  result = (8-8*J-7*J^2+2*J^3+J^4)/(2*sp^5)-2*2/sp^5;
  Re[N[result, 600] ]
];

n4ABBA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = (sp (12 mA^4+6 (-4+J+J^2) mA^2 sp+(12-8 J-7 J^2+2 J^3+J^4) sp^2))/(4 (mA^2-sp)^8);
  Re[N[result, 600] ]
];

n4BBAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = -(((-(sp/(4 mA^2-sp)))^(1/4) LegendreP[J,2,2,(-2 mA^2+sp)/Sqrt[sp (-4 mA^2+sp)]])/(2 mA^4 (mA^2-sp)^3));
  Re[N[result, 600] ]
];

n4ABAB[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = (sp (12 LegendreP[J,(-mA^4-2 mA^2 sp+sp^2)/(mA^2-sp)^2]+(6 (-mA^2+sp) LegendreP[J,1,2,(-mA^4-2 mA^2 sp+sp^2)/(mA^2-sp)^2])/(mA^2 Sqrt[1-(2 mA^2)/sp])-((mA^2-sp)^2 sp LegendreP[J,2,2,(-mA^4-2 mA^2 sp+sp^2)/(mA^2-sp)^2])/(2 mA^6-mA^4 sp)))/(4 (mA^2-sp)^5 (-mA^2+sp));
  Re[N[result, 600] ]
];


X52ABBA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA, result}, 
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = (sp (5 mA^2+(-5+2 J+2 J^2) sp))/(mA^2-sp)^8;
  Re[N[result, 600] ]
];

X52BBAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA,result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = ((Sqrt[sp/(4 mA^2-sp)] LegendreP[J,3,2,(2 mA^2-sp)/Sqrt[sp (-4 mA^2+sp)]])/(48 mA^6 (mA^2-sp)^3));
  Re[N[result, 600] ]
];

X52ABAB[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA,result},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  result = 1/(3 (mA^2-sp)^6 (-mA^2+sp)) sp (15 LegendreP[J,(2 (2 mA^2-sp) sp)/(mA^2-sp)^2]+(18 (mA^2-sp) sp LegendreP[J,1,2,(2 (2 mA^2-sp) sp)/(mA^2-sp)^2])/Sqrt[mA^8-4 mA^6 sp-10 mA^4 sp^2+12 mA^2 sp^3-3 sp^4]+(9 (mA^2-sp)^2 sp^2 LegendreP[J,2,2,(2 (2 mA^2-sp) sp)/(mA^2-sp)^2])/(mA^8-4 mA^6 sp-10 mA^4 sp^2+12 mA^2 sp^3-3 sp^4)+(2 (mA^2-sp)^3 sp^3 LegendreP[J,3,2,(2 (2 mA^2-sp) sp)/(mA^2-sp)^2])/(mA^8-4 mA^6 sp-10 mA^4 sp^2+12 mA^2 sp^3-3 sp^4)^(3/2));
  Re[N[result, 600] ]
];


X52AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

X52BBBB[x_?NumericQ, J_?IntegerQ] := Module[{sp},
  sp = N[1/(1-x), 600];
  N[-((J (1+J) (150+J (1+J) (-43+2 J (1+J))))/(36 sp^6)), 600]
];

X62AAAA[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

X62BBBB[x_?NumericQ, J_?IntegerQ] := Module[{sp},
  sp = N[1/(1-x), 600];
  N[-(((-3+J) J (1+J) (4+J) (204+J (1+J) (-32+J+J^2)))/(288 sp^7)), 600]
];

(* Large J limit k = 6, q = 2 *)
LargeJAAAA[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];    (* FIX 1: exact rational *)
  N[-(((1/(-4 mA^2 + sp))^(15/2) (-32 mA^6 + 24 mA^4 sp - 6 mA^2 sp^2 + sp^3))/(288 sp^(5/2))), 600]
];

LargeJBBBB[x_?NumericQ] := Module[{sp},
  sp = N[1/(1-x), 600];
  N[-(1/(288 sp^7)), 600]
];

Jmax = 60;
Jlist = Range[0, Jmax, 2];

M0[x_?NumericQ,J_?IntegerQ] := {
	{g20[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M1[x_?NumericQ,J_?IntegerQ] := {
	{g31[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M2[x_?NumericQ,J_?IntegerQ] := {
  {g4[x,J],0,0},
  {0,0,0},
  {0,0,0}
};

M41[x_?NumericQ, J_?IntegerQ] :={
	{n4AAAA[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M42[x_?NumericQ, J_?IntegerQ] :={
	{0,0,0},
	{0,n4BBBB[x,J],0},
	{0,0,0}
};

M43[x_?NumericQ, J_?IntegerQ] :={
  {0,1/2*n4BBAA[x,J],0},
  {1/2*n4BBAA[x,J],0,0},
  {0,0,n4ABAB[x,J]+n4ABBA[x,J]}
};

M51[x_?NumericQ, J_?IntegerQ] :={
	{X52AAAA[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M52[x_?NumericQ, J_?IntegerQ] :={
  {0,0,0},
	{0,X52BBBB[x,J],0},
	{0,0,0}
};

M53[x_?NumericQ, J_?IntegerQ] :={
  {0,1/2*X52BBAA[x,J],0},
  {1/2*X52BBAA[x,J],0,0},
  {0,0,X52ABAB[x,J]+X52ABBA[x,J]}
};

M61[x_?NumericQ, J_?IntegerQ] :={
	{X62AAAA[x,J],0,0},
	{0,0,0},
	{0,0,0}
};

M62[x_?NumericQ, J_?IntegerQ] :={
  {0,0,0},
	{0,X62BBBB[x,J],0},
	{0,0,0}
};

M6j1[x_?NumericQ] :={
	{LargeJAAAA[x],0,0},
	{0,0,0},
	{0,0,0}
}

M6j2[x_?NumericQ] :={
	{0,0,0},
	{0,LargeJBBBB[x],0},
	{0,0,0}
}

(* null *)
N0[x_?NumericQ] :={
	{0,0,0},
	{0,0,0},
	{0,0,0}
};

f11List ={
	Function[{x,J}, M0[x,J][[1,1]]],
	Function[{x,J}, M1[x,J][[1,1]]],
  Function[{x,J}, M2[x,J][[1,1]]],
	Function[{x,J}, M41[x,J][[1,1]]],
  Function[{x,J}, M42[x,J][[1,1]]],
  Function[{x,J}, M43[x,J][[1,1]]],
	Function[{x,J}, M51[x,J][[1,1]]],
  Function[{x,J}, M52[x,J][[1,1]]],
  Function[{x,J}, M53[x,J][[1,1]]],
  Function[{x,J}, M61[x,J][[1,1]]],
  Function[{x,J}, M62[x,J][[1,1]]]
};

f22List ={
	Function[{x,J}, M0[x,J][[2,2]]],
	Function[{x,J}, M1[x,J][[2,2]]],
  Function[{x,J}, M2[x,J][[2,2]]],
	Function[{x,J}, M41[x,J][[2,2]]],
  Function[{x,J}, M42[x,J][[2,2]]],
  Function[{x,J}, M43[x,J][[2,2]]],
	Function[{x,J}, M51[x,J][[2,2]]],
  Function[{x,J}, M52[x,J][[2,2]]],
  Function[{x,J}, M53[x,J][[2,2]]],
  Function[{x,J}, M61[x,J][[2,2]]],
  Function[{x,J}, M62[x,J][[2,2]]]
};

f33List = {
	Function[{x,J}, M0[x,J][[3,3]]],
	Function[{x,J}, M1[x,J][[3,3]]],
  Function[{x,J}, M2[x,J][[3,3]]],
	Function[{x,J}, M41[x,J][[3,3]]],
  Function[{x,J}, M42[x,J][[3,3]]],
  Function[{x,J}, M43[x,J][[3,3]]],
	Function[{x,J}, M51[x,J][[3,3]]],
  Function[{x,J}, M52[x,J][[3,3]]],
  Function[{x,J}, M53[x,J][[3,3]]],
  Function[{x,J}, M61[x,J][[3,3]]],
  Function[{x,J}, M62[x,J][[3,3]]]
};

f12List ={
	Function[{x,J}, M0[x,J][[1,2]]],
	Function[{x,J}, M1[x,J][[1,2]]],
  Function[{x,J}, M2[x,J][[1,2]]],
	Function[{x,J}, M41[x,J][[1,2]]],
  Function[{x,J}, M42[x,J][[1,2]]],
  Function[{x,J}, M43[x,J][[1,2]]],
	Function[{x,J}, M51[x,J][[1,2]]],
  Function[{x,J}, M52[x,J][[1,2]]],
  Function[{x,J}, M53[x,J][[1,2]]],
  Function[{x,J}, M61[x,J][[1,2]]],
  Function[{x,J}, M62[x,J][[1,2]]]
};

f21List = f12List;

f13List = {
	Function[{x,J}, M0[x,J][[1,3]]],
	Function[{x,J}, M1[x,J][[1,3]]],
  Function[{x,J}, M2[x,J][[1,3]]],
	Function[{x,J}, M41[x,J][[1,3]]],
  Function[{x,J}, M42[x,J][[1,3]]],
  Function[{x,J}, M43[x,J][[1,3]]],
	Function[{x,J}, M51[x,J][[1,3]]],
  Function[{x,J}, M52[x,J][[1,3]]],
  Function[{x,J}, M53[x,J][[1,3]]],
  Function[{x,J}, M61[x,J][[1,3]]],
  Function[{x,J}, M62[x,J][[1,3]]]
};
f31List = f13List;

f23List = {
	Function[{x,J}, M0[x,J][[2,3]]],
	Function[{x,J}, M1[x,J][[2,3]]],
  Function[{x,J}, M2[x,J][[2,3]]],
	Function[{x,J}, M41[x,J][[2,3]]],
  Function[{x,J}, M42[x,J][[2,3]]],
  Function[{x,J}, M43[x,J][[2,3]]],
	Function[{x,J}, M51[x,J][[2,3]]],
  Function[{x,J}, M52[x,J][[2,3]]],
  Function[{x,J}, M53[x,J][[2,3]]],
  Function[{x,J}, M61[x,J][[2,3]]],
  Function[{x,J}, M62[x,J][[2,3]]]
};
f32List = f23List;

j11List = {
	Function[{x}, N0[x][[1,1]]],
	Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
	Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, N0[x][[1,1]]],
  Function[{x}, M6j1[x][[1,1]]],
	Function[{x}, M6j2[x][[1,1]]]
};
j22List = {
	Function[{x}, N0[x][[2,2]]],
	Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
	Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, N0[x][[2,2]]],
  Function[{x}, M6j1[x][[2,2]]],
	Function[{x}, M6j2[x][[2,2]]]
};
j33List = {
  Function[{x}, N0[x][[3,3]]],
	Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
	Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, N0[x][[3,3]]],
  Function[{x}, M6j1[x][[3,3]]],
	Function[{x}, M6j2[x][[3,3]]]
};

j12List = {
  Function[{x}, N0[x][[1,2]]],
	Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
	Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, N0[x][[1,2]]],
  Function[{x}, M6j1[x][[1,2]]],
	Function[{x}, M6j2[x][[1,2]]]
};
j21List = j12List;

j13List = {
  Function[{x}, N0[x][[1,3]]],
	Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
	Function[{x}, N0[x][[1,3]]],
	Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
  Function[{x}, N0[x][[1,3]]],
  Function[{x}, M6j1[x][[1,3]]],
	Function[{x}, M6j2[x][[1,3]]]
};
j31List = j13List;

j23List = {
  Function[{x}, N0[x][[2,3]]],
	Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
	Function[{x}, N0[x][[2,3]]],
	Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
  Function[{x}, N0[x][[2,3]]],
  Function[{x}, M6j1[x][[2,3]]],
	Function[{x}, M6j2[x][[2,3]]]
};
j32List = j23List;


(* large J limit: 0& for g20,g31,n4,X52; LargeJ for X53 *)
(* extraTriplet = {0&, 0&, 0&, LargeJ}; *)

(* optimal upper bound *)
norm = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
obj  = {-1, 5.30625916733434, 0, 0, 0, 0, 0, 0, 0, 0, 0};


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

  Print["  Regular blocks : ", Length[samplePoints] * Length[Jlist]];
  Print["  Extra blocks   : ", Length[polsExtra]];
  Print["  Total blocks   : ", Length[samplePoints] * Length[Jlist] + Length[polsExtra]];

  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Join[Flatten[polsRegular], polsExtra]],
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

    Print["=== test9.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    Null
  ]
];
