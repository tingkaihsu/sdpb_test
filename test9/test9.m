(* pmp generator for AAAA scattering *)
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
    (* functionValues has structure {{{ {f1(xᵢ)}, {f2(xᵢ)} }}} — 3 structural
       wrapper levels, with depth-3 entries being the degree-0 coefficient lists.
       toJsonNestedNumberArray stringifies all numeric leaves in-place. *)
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
(* let the mass be m = 0.2 so that 4m^2 < M^2 = 1 where M  = 1 to infinity *)

maVal = SetPrecision[0.040, 50];

(* dispersion representation of Wilson coefficients *)
(* All functions now precompute sp and mA numerically with N[...,50] *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[1/2*Sqrt[sp/(sp - 4*mA^2)] * (sp^(-3) + (-4*mA^2 + sp)^(-3)), 50]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(-Sqrt[sp/(sp - 4*mA^2)] * ((-3 + J*(1 + J)*(-4*mA^2 + sp)^3*(sp^(-3) + (-4*mA^2 + sp)^(-3)))/(-4*mA^2 + sp)^4)), 50]
];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(4*mA - sp)^(-5) + (4 - (-2 + J)*J*(1 + J)*(3 + J))/(4*sp^5) + (2*J*(1 + J))/(sp*(-4*mA + sp)^4) - ((-1 + J)*J*(1 + J)*(2 + J))/(4*sp^2*(-4*mA + sp)^3), 50]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 50]
];

X62[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(-2304*mA^4 + 1152*mA^2*sp + (-72 + J*(1 + J)*(-18 + J + J^2))*sp^2))/sp^5))/(576*(-4*mA^2 + sp)^7), 50]
];

X72[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(-47 + J + J^2))/(-4*mA^2 + sp)^8) + ((-1 + J)*(2 + J)*(((-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*sp^3)/(4*mA^2 - sp)^5 + 3600/(-4*mA^2 + sp)^2))/sp^6))/14400, 50]
];

X82[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(((-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(-38 + J + J^2))/(4*mA^2 - sp)^9 + ((-1 + J)*(2 + J)*(-(((-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*sp^4)/(-4*mA^2 + sp)^6) + 129600/(-4*mA^2 + sp)^2))/sp^7))/518400, 50]
];

X92[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(2754 + J*(1 + J)*(-119 + J + J^2)))/(-4*mA^2 + sp)^10) + ((-1 + J)*(2 + J)*(((-6 + J)*(-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*(7 + J)*sp^5)/(4*mA^2 - sp)^7 + 6350400/(-4*mA^2 + sp)^2))/sp^8))/25401600, 50]
];

X102[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-(((-10 + J)*(-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(11 + J)*(2232 + (-10 + J)*J*(1 + J)*(11 + J)))/(-4*mA^2 + sp)^11) + ((-1 + J)*(2 + J)*(-(((-7 + J)*(-6 + J)*(-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*(7 + J)*(8 + J)*sp^6)/(-4*mA^2 + sp)^8) + 406425600/(-4*mA^2 + sp)^2))/sp^9))/1625702400, 50]
];

(* Large J limit *)
LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 50];
  mA = N[maVal, 50];
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-(-4*mA^2 + sp)^(-11) - 1/(sp^3*(-4*mA^2 + sp)^8)))/1625702400, 50]
];

Jmax = 40;
Jlist = Range[0, Jmax, 2];

fList = {g20, g31, n4, X52, X62, X72, X82, X92, X102};

(* large J limit *)
(* 0& is a constant function of 0 *)

extraTriplet = {0&, 0&, 0&, 0&, 0&, 0&, 0&, 0&, LargeJ};

norm = {0, 1, 0, 0, 0, 0, 0, 0, 0};
obj = {-1, 0, 0, 0, 0, 0, 0, 0, 0};

testNumericalSDP[spFile_String, jsonFile_String, prec_:200] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, polsRegular},

  (* --- Read and parse sampling_points.txt --- *)
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
  Print["  extraTriplet: ", extraTriplet, "  (J\[Rule]\[Infinity] limit)"];

  (* Scalings = prefactor DampedRational[1,{},1/E,x] evaluated at xi = e^{-xi} *)
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  (* --- Regular blocks: one per (xi, Jj) pair.
     Polynomials nesting:  {{ Table[{fk(xi,Jj)}, {k,3}] }}
       {{ ... }}  ← JSON levels 1 and 2  (column list / row, each size 1)
       Table[...] ← JSON level 3: polynomial vector, one entry per fList[[k]]
       {fk(xi,Jj)} ← JSON level 4: degree-0 coefficient list (1 element) --- *)
  polsRegular = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      (* TWO outer braces = JSON levels 1 & 2.
         Table produces level-3 list: {{f1(xi,Jj)}, {f2(xi,Jj)}, {f3(xi,Jj)}}.
         Each {fList[[k]][...]} is the level-4 coefficient list (degree-0: 1 entry). *)
      "polynomials" -> {{ Table[
        {SetPrecision[fList[[k]][samplePoints[[i]], Jlist[[j]]], prec]},
        {k, Length[fList]}
      ] }}
    |>],
    {i, Length[samplePoints]}, {j, Length[Jlist]}
  ];

  (* large J limit *)
  polsExtra = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[extraTriplet[[k]][samplePoints[[i]]], prec]},
        {k, Length[extraTriplet]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  Print["  Regular blocks : ", Length[samplePoints] * Length[Jlist]];
  Print["  Extra blocks   : ", Length[polsExtra]];
  Print["  Total blocks   : ", Length[samplePoints] * Length[Jlist] + Length[polsExtra]];

  (* Flatten polsRegular (2D Table → flat list) *)
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
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 200];

    Print["=== text9.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    (* Loaded with << as a library — do nothing. *)
    Null
  ]
];
