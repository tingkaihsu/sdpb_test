(* New head for purely numerical data *)
ClearAll[NumericalPositiveMatrixWithPrefactor];

(* Helper: recursively stringify all numeric leaves *)
toJsonNestedNumberArray[expr_, prec_] := expr /. n_?NumericQ :> toJsonNumber[n, prec];

toJsonObject[NumericalPositiveMatrixWithPrefactor[pmp_?AssociationQ], prec_, getSampleDataFn_:Function[<||>]] :=
Module[
  {sampleData, samplePoints, sampleScalings, functionValues, basisValues,
   prefactor, reducedPrefactor},

  sampleData = getSampleDataFn[NumericalPositiveMatrixWithPrefactor[pmp], prec];

  prefactor = Lookup[pmp, "prefactor", 1];
  reducedPrefactor = Lookup[pmp, "reducedPrefactor", prefactor];

  samplePoints   = Lookup[sampleData, "samplePoints",   Lookup[pmp, "samplePoints",   Missing[]]];
  sampleScalings  = Lookup[sampleData, "sampleScalings",  Lookup[pmp, "sampleScalings",  Missing[]]];
  functionValues = Lookup[pmp, "polynomials", Missing[]];
  basisValues    = Lookup[sampleData, "basisValues", Lookup[pmp, "basisValues", Missing[]]];

  DeleteMissing @ <|
    "prefactor" -> toJsonDampedRational[prefactor, prec],
    "reducedPrefactor" -> toJsonDampedRational[reducedPrefactor, prec],
    "samplePoints" -> toJsonNumberArray[samplePoints, prec],
    "sampleScalings" -> toJsonNumberArray[sampleScalings, prec],
    "polynomials" -> If[MissingQ[functionValues], Missing[], toJsonNestedNumberArray[functionValues, prec]],
    "basisValues" -> If[MissingQ[basisValues], Missing[], toJsonNestedNumberArray[basisValues, prec]]
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
    "objective" -> toJsonNumberArray[objective, prec],
    "normalization" -> toJsonNumberArray[normalization, prec],
    "PositiveMatrixWithPrefactorArray" ->
      Table[toJsonObject[pmp, prec, getSampleDataFn], {pmp, positiveMatricesWithPrefactors}]
  |>
];


<< "../SDPB.m";

(* Example numeric functions; in a real application these can be NIntegrate-based
   or any black-box numerical evaluator. *)
f1[J_?NumericQ, x_?NumericQ] := 1 + x^4;
f2[J_?NumericQ, x_?NumericQ] := x^4/12 + x^2 + J;

myNumericSampleData[NumericalPositiveMatrixWithPrefactor[pmp_], prec_] := <|
  "samplePoints" -> SetPrecision[{0.1, 0.5, 1.2, 2.8, 6.0}, prec],
  "sampleScalings" -> SetPrecision[{0.2, 0.6, 1.1, 1.7, 2.4}, prec]
|>;

testNumericalSDP[jsonFile_, prec_:200] := Module[
  {Jlist, pols, norm, obj, samplePoints},

  Jlist = Range[0, 40, 2];
  samplePoints = SetPrecision[{0.1, 0.5, 1.2, 2.8, 6.0}, prec];

  pols = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor" -> DampedRational[1, {}, 1/E, x],
      "samplePoints" -> samplePoints,
      "sampleScalings" -> SetPrecision[{0.2, 0.6, 1.1, 1.7, 2.4}, prec],
      "polynomials" -> Table[
        {
          f1[J, xi],
          f2[J, xi]
        },
        {xi, samplePoints}
      ]
    |>],
    {J, Jlist}
  ];

  norm = {1, 0};
  obj  = {0, -1};

  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, pols],
    prec,
    myNumericSampleData
  ]
];

testNumericalSDP["numeric_pmp.json", 200];