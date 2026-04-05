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
f1[x_?NumericQ] := 1 + x^4;
f2[x_?NumericQ] := x^4/12 + x^2;

myNumericSampleData[NumericalPositiveMatrixWithPrefactor[pmp_], prec_] := <|
  "samplePoints" -> SetPrecision[{0.06, 0.57, 1.63, 3.42, 6.49}, prec],
  "sampleScalings" -> SetPrecision[{0.94, 0.57, 0.20, 0.03, 0.001}, prec]
|>;

testNumericalSDP[jsonFile_, prec_:200] := Module[
  {pols, norm, obj, samplePoints},

  samplePoints = SetPrecision[{0.06, 0.57, 1.63, 3.42, 6.49}, prec];

  pols =
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor" -> DampedRational[1, {}, 1/E, x],
      "samplePoints" -> samplePoints,
      "sampleScalings" -> SetPrecision[{0.94, 0.57, 0.20, 0.03, 0.001}, prec],
      "polynomials" -> {
        {
          {
            Table[f1[xi], {xi, samplePoints}],
            Table[f2[xi], {xi, samplePoints}]
          }
        }
      }
    |>];

  norm = {1, 0};
  obj  = {0, -1};

  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, {pols}],
    prec,
    myNumericSampleData
  ]
];

testNumericalSDP["numeric_pmp.json", 2];