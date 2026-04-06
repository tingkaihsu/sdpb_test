(* New head for purely numerical data *)
ClearAll[NumericalPositiveMatrixWithPrefactor];

(* Helper: recursively stringify all numeric leaves *)
toJsonNestedNumberArray[expr_, prec_] := expr /. n_?NumericQ :> toJsonNumber[n, prec];

(*
  SDPB manual (Section 3.1) is definitive:
      polynomial c0 + c1*x + ... + cd*x^d  <=>  ["c0", ..., "cd"]
  The `polynomials` field ALWAYS requires polynomial coefficient lists.
  There is no "sample-values mode."

  The analytic path (SDPB.m) has symbolic polynomial expressions in pmp[["polynomials"]]
  and extracts coefficients directly via safeCoefficientList[#, x].

  The numerical path (here) stores function VALUES evaluated at sample points instead,
  because the functions are black-box evaluators (e.g. NIntegrate-based) with no
  accessible polynomial form. To produce coefficient lists we must reconstruct the
  polynomial via Lagrange interpolation.

  Previous attempt used the user-specified `prec` as the working precision for the
  interpolation computation. This fails because the 5 sample points over [0.06, 6.5]
  give a Vandermonde matrix with condition number ~10^8: one needs prec + 8 digits of
  working precision minimum just to recover 2-digit coefficients. The fix: use
  Max[prec + 50, 100] digits internally, then round to `prec` at output. This is
  stable for all practical values of prec >= 2.
*)

interpolatedCoefficients[samplePoints_List, vals_List, prec_] :=
  Module[
    {workPrec = Max[prec + 50, 100], poly},
    (* Lift both sample points and values to workPrec before interpolating.
       This ensures the divided-difference computation does not lose information
       to rounding, regardless of how low `prec` is. *)
    poly = Expand @ InterpolatingPolynomial[
      Transpose[{
        SetPrecision[samplePoints, workPrec],
        SetPrecision[vals,         workPrec]
      }],
      x
    ];
    (* Extract coefficient list, then round to the user-requested precision. *)
    SetPrecision[safeCoefficientList[poly, x], prec]
  ];

toJsonObject[NumericalPositiveMatrixWithPrefactor[pmp_?AssociationQ], prec_, getSampleDataFn_:Function[<||>]] :=
Module[
  {sampleData, samplePoints, sampleScalings, functionValues, basisValues,
   prefactor, reducedPrefactor},

  sampleData = getSampleDataFn[NumericalPositiveMatrixWithPrefactor[pmp], prec];

  prefactor        = Lookup[pmp, "prefactor",        1];
  reducedPrefactor = Lookup[pmp, "reducedPrefactor", prefactor];

  (* samplePoints are taken from sampleData first (the authoritative source used
     for evaluating the functions), falling back to the pmp association. *)
  samplePoints   = Lookup[sampleData, "samplePoints",  Lookup[pmp, "samplePoints",  Missing[]]];
  sampleScalings = Lookup[sampleData, "sampleScalings", Lookup[pmp, "sampleScalings", Missing[]]];
  functionValues = Lookup[pmp, "polynomials", Missing[]];
  basisValues    = Lookup[sampleData, "basisValues", Lookup[pmp, "basisValues", Missing[]]];

  DeleteMissing @ <|
    "prefactor"        -> toJsonDampedRational[prefactor,        prec],
    "reducedPrefactor" -> toJsonDampedRational[reducedPrefactor, prec],
    "samplePoints"     -> toJsonNumberArray[samplePoints,   prec],
    "sampleScalings"   -> toJsonNumberArray[sampleScalings, prec],

    (* Reconstruct polynomial coefficient lists from sampled values via Lagrange
       interpolation. Map at depth {3} mirrors the analytic path in SDPB.m:
         Map[toJsonNumberArray[safeCoefficientList[#,x],prec]&, pmp[["polynomials"]], {3}]
       In both cases the level-3 elements are the per-entry polynomial data:
         analytic: a symbolic expression like  1 + x^4
         numeric:  a list of sampled values    {f(x0), f(x1), ..., f(xd)} *)
    "polynomials" -> If[MissingQ[functionValues] || MissingQ[samplePoints],
      Missing[],
      Map[
        Function[vals,
          toJsonNumberArray[
            interpolatedCoefficients[samplePoints, vals, prec],
            prec
          ]
        ],
        functionValues,
        {3}
      ]
    ],

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
    "objective"     -> toJsonNumberArray[objective,     prec],
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
  "samplePoints"   -> SetPrecision[{0.06, 0.57, 1.63, 3.42, 6.49}, prec],
  "sampleScalings" -> SetPrecision[{0.94, 0.57, 0.20, 0.03, 0.001}, prec]
|>;

testNumericalSDP[jsonFile_, prec_:200] := Module[
  {pols, norm, obj, samplePoints},

  samplePoints = SetPrecision[{0.06, 0.57, 1.63, 3.42, 6.49}, prec];

  pols =
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> samplePoints,
      "sampleScalings" -> SetPrecision[{0.94, 0.57, 0.20, 0.03, 0.001}, prec],
      (* Store function values at sample points.
         toJsonObject will reconstruct polynomial coefficient lists via interpolation. *)
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

testNumericalSDP["numeric_pmp.json"];