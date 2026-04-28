(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* AAAA scattering *)
Mlow[s_, t_, Nmax_Integer] :=
    gAAB^2 * (1/s + 1/t + 1/u) +
    gAAA^2 * (1/(s-1) + 1/(t-1) + 1/(u-1)) +
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s)^ab[[1]]
            * (t)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4 - s - t};


(* dispersion representation of Wilson coefficients *)

v[l_, q_] := Product[l*(l + 1) - a*(a - 1), {a, 1, q}] / (Factorial[q])^2;

(* Print["v[l,0] = ", v[l, 0]//FullSimplify]; *)
(* 1 *)

(* Print["v[l,1] = ", v[l, 1]//FullSimplify]; *)
(* l(l+1) *)

(* Null constraints *)

kernelDirect[k_Integer, q_Integer, u_, sp_, mA_] :=
  1/u (1/(u^q sp^(k - q + 1)) - 1/(u^(k - q) sp^(q + 1)));

kernelCross[k_Integer, q_Integer, u_, sp_, mA_] :=
  1/u (1/(u^q (4mA^2 -sp -u)^(k - q + 1)) - 1/(u^(k - q) (4mA^2 -sp -u)^(q + 1)));



directPiece[k_Integer, q_Integer, J_, t_, sp_, mA_] :=
  kernelDirect[k, q, t, sp, mA] * sp * LegendreP[J, 1 + 2 t/(sp-4mA^2)];

crossPiece[k_Integer, q_Integer, J_, t_, sp_, mA_] :=
  kernelCross[k, q, t, sp, mA] * sp * LegendreP[J, 1 + 2 t/(sp-4mA^2)];

directRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[directPiece[8, 2, J, t, sp, mA], {t, 0, -1}]
];

crossRes = Assuming[J ∈ Integers && J >= 0,
  FullSimplify @ SeriesCoefficient[crossPiece[8, 2, J, t, sp, mA], {t, 0, -1}]
];

combined = FullSimplify[directRes - crossRes];

(* Print["Combined Piece: ", combined]; *)

(* k, q *)
(* 5, 2 *)
(* (J*(1 + J)*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)*sp^4) - (-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp)))/(36*sp^3*(-4*mA^2 + sp)^6) *)

(* 6, 2 *)
(* (J*(1 + J)*(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*sp^5)/(4*mA^2 - sp)^7 + ((-2 + J + J^2)*(2304*mA^4 - 1152*mA^2*sp + (72 - J*(1 + J)*(-18 + J + J^2))*sp^2))/(-4*mA^2 + sp)^4))/(576*sp^4) *)

(* 7, 2 *)
(* (J*(1 + J)*(-(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(-47 + J + J^2)*sp^6)/(-4*mA^2 + sp)^8) + (-1 + J)*(2 + J)*(((-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*sp^3)/(4*mA^2 - sp)^5 + 3600/(-4*mA^2 + sp)^2)))/(14400*sp^5) *)

(* 8, 2 *)
(* (J*(1 + J)*(((-8 + J)*(-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(9 + J)*(-38 + J + J^2)*sp^7)/(4*mA^2 - sp)^9 + (-1 + J)*(2 + J)*(-(((-5 + J)*(-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*(6 + J)*sp^4)/(-4*mA^2 + sp)^6) + 129600/(-4*mA^2 + sp)^2)))/(518400*sp^6) *)

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

ma = 0.2;

f1[x_?NumericQ, J_?IntegerQ] := Sqrt[ M^2/ (M^2-4*mA^2) ] * v[J, 0] * (1/M^2)^2/.{M -> Sqrt[1/(1-x)], mA -> ma};

f2[x_?NumericQ, J_?IntegerQ] := Sqrt[ M^2/ (M^2-4*mA^2) ] * v[J, 1] * (1/M^2)^2 * 1/(M^2-4*mA^2)/.{M -> Sqrt[1/(1-x)], mA -> ma};

X52[x_?NumericQ, J_?IntegerQ] := (J*(1 + J)*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)*sp^4) - (-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp)))/(36*sp^3*(-4*mA^2 + sp)^6)/.{mA -> ma}/.{sp -> 1/(1-x)};
X62[x_?NumericQ, J_?IntegerQ] := (J*(1 + J)*(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*sp^5)/(4*mA^2 - sp)^7 + ((-2 + J + J^2)*(2304*mA^4 - 1152*mA^2*sp + (72 - J*(1 + J)*(-18 + J + J^2))*sp^2))/(-4*mA^2 + sp)^4))/(576*sp^4)/.{mA -> ma}/.{sp -> 1/(1-x)};
X72[x_?NumericQ, J_?IntegerQ] := (J*(1 + J)*(-(((-6 + J)*(-4 + J)*(-2 + J)*(3 + J)*(5 + J)*(7 + J)*(-47 + J + J^2)*sp^6)/(-4*mA^2 + sp)^8) + (-1 + J)*(2 + J)*(((-4 + J)*(-3 + J)*(-2 + J)*(3 + J)*(4 + J)*(5 + J)*sp^3)/(4*mA^2 - sp)^5 + 3600/(-4*mA^2 + sp)^2)))/(14400*sp^5)/.{mA -> ma}/.{sp -> 1/(1-x)};

Jmax = 40;
Jlist = Range[0, Jmax, 2];

fList = {f1, f2, X52, X62, X72};

norm = {0, 1, 0, 0, 0};
obj = {-1, 0, 0, 0, 0};

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

  Print["  Regular blocks : ", Length[samplePoints] * Length[Jlist]];
  Print["  Total blocks   : ", Length[samplePoints] * Length[Jlist]];

  (* Flatten polsRegular (2D Table → flat list) *)
  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Flatten[polsRegular]],
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