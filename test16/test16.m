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

mA = N[1/1000, 600];

Print["mA = ", mA]

(* g2 and g3 in AB -> AB scattering *)
cList[s_,J_] := {
    s/(mA^2-s)^4,
    (s (3 mA^2+(-3+2 J (1+J)) s))/(2 (mA^2-s)^6)
};

nullList = {-1, -1, -1, -1};

list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];

Nlist[n_,s_,mA_,J_] := 

Jmax = 70;
Jlist = Range[0, Jmax, 2];


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