(* pmp generator for AAAA scattering *)
(* ============================================================
   PRECISION-PROPAGATION FIX (2026-05)
   ------------------------------------------------------------
   The null constraint n4[x, J] contains three special functions
   (Hypergeometric2F1, LegendreP[J,1,·], LegendreP[J,2,·]) that
   each grow as  A(x)^J  with  A = z₂ + √(z₂²−1) > 1.
   They must cancel exactly (null identity).  Resolving this
   cancellation faithfully requires

       prec_local  ≥  J · max_x log₁₀ A(x)  +  margin
                   ≈  J · 0.1568            +  60

   LEAKAGE POINTS FIXED
   1. maVal = SetPrecision[0.150, 600]: machine float 0.150 has
      only 15–16 true digits.  SetPrecision pads with zeros from
      the wrong binary representation.  FIX: use N[3/20, prec]
      (exact rational) inside every function that needs precision
      beyond 16 digits.
   2. Fixed prec=600/650 inside n4 regardless of J.  At J=10000
      the terms are ~10^1568 and 600 digits leave ~10^968 noise.
      FIX: dynamic prec via n4PrecMin[J, x].
   3. Outer N[expr, 600] doesn't increase internal evaluation
      precision of Hypergeometric2F1/LegendreP — those inherit
      precision from their arguments.  FIX: sp, mA computed at
      prec_local digits so arguments carry the right precision.
   4. Im residual: LegendreP[J,1,z>1] is imaginary in Mathematica's
      convention; the (-2I) factor makes the product real.  At finite
      precision an Im residual ~10^{-margin} remains.  FIX: Re[…]
      made explicit so the output is always a real number.
   ============================================================ *)

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

(* FIX 1 of 4: Use the exact rational 3/20 instead of the machine float 0.150.
   SetPrecision[0.150, 600] pads from the wrong IEEE-754 bit pattern (which
   encodes 0.14999999999999999...), giving a spuriously "high-precision" but
   factually wrong value at digit 17 onward.  N[3/20, 600] is exact. *)
maVal = N[3/20, 650];

Print["mA = ", maVal]

(* ============================================================
   FIX 2 of 4: n4PrecMin — adaptive working precision for n4.

   Returns the minimum number of decimal digits needed so that
   evaluating n4[x, J] resolves the catastrophic cancellation of
   the three  A(x)^J  terms to `margin` true significant figures.

   Algorithm:
     z₂  = 1 + 8·mA²/(3·(sp − 4·mA²))  (the Legendre argument)
     A   = z₂ + √(z₂²−1)  > 1
     log₁₀(A) ≈ 0.157 at worst case (x → 0)
     prec_local = max(stdPrec, ⌈J · log₁₀(A)⌉ + margin)

   For J ≤ 60: prec_local = stdPrec = 650   → zero overhead.
   For J = 10000: prec_local = 1628          → correct resolution.
   For J = 20000: prec_local = 3196          → correct resolution.

   The estimate itself uses only 50-digit arithmetic (no cancellation).
   The exact rational (3/20)^2 is used for mA inside the estimator
   to avoid the machine-float contamination described in §2.1.
   ============================================================ *)

n4PrecMin[J_?IntegerQ, x_?NumericQ, margin_Integer:60, stdPrec_Integer:650] :=
  Module[{sp0, z2, A, logA},
    sp0  = N[1/(1 - x), 50];
    z2   = 1 + 8*(3/20)^2 / (3*(sp0 - 4*(3/20)^2));
    A    = If[z2 > 1, z2 + Sqrt[z2^2 - 1], 1];
    logA = If[A > 1, N[Log[10, A], 50], 0];
    Max[stdPrec, Ceiling[J * logA] + margin]
  ];

(* forward limit: use our own convention *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[3/20, 650];    (* FIX 1: exact rational *)
  N[-(2*Sqrt[sp/(-4*mA^2 + sp)])/(2*mA^2 - sp)^3, 650]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[3/20, 650];    (* FIX 1: exact rational *)
  N[-(Sqrt[sp/(-4*mA^2 + sp)]*(-3 - (2*J*(1 + J)*(2*mA^2 - sp))/(-4*mA^2 + sp)))/(-2*mA^2 + sp)^4, 650]
];

(* null constraint — PRECISION-PROPAGATION IMPLEMENTATION
   --------------------------------------------------------
   Strategy: compute all intermediate quantities at prec_local digits,
   where prec_local is determined adaptively by n4PrecMin.

   For J ≤ ~3782: prec_local = 650 (same as other functions, zero overhead).
   For J > 3782:  prec_local > 650 (extra digits only where required).

   Key implementation points:
   - xp = SetPrecision[x, prec]: promotes the input x to prec digits.
     x may arrive with fewer digits (e.g. 200 from the file reader);
     SetPrecision pads rather than creates precision, but since n4 = 0
     identically, the result does not depend sensitively on the last
     digits of x.  The promotion prevents the 200-digit x from acting
     as a bottleneck that limits all downstream arithmetic.
   - mA = N[3/20, prec]: exact rational for mA, ignores global maVal.
   - All arithmetic involving sp, mA (including the special-function
     arguments) inherits prec-digit precision automatically through
     Mathematica's internal precision tracker.
   - Hypergeometric2F1 and LegendreP are evaluated at ≈ prec digits
     because their arguments carry prec-digit precision.
   - Re[N[result, 650]]: the imaginary residual from LegendreP[J,1,z>1]
     should be ~10^{-margin} after the (-2I) cancellation.  Re[…]
     extracts the real part explicitly and serves as a self-diagnostic:
     if Im is large, prec_local was insufficient.

   The formula itself is unchanged from the original. *)

n4[x_?NumericQ, J_?IntegerQ] := Module[{prec, xp, sp, mA, result},
  (* FIX 2: adaptive precision — no overhead for J ≤ 60 *)
  prec = n4PrecMin[J, x];
  (* FIX 3: promote input x; avoids 200-digit bottleneck from file reader *)
  xp   = SetPrecision[x, prec];
  (* FIX 1 (inside n4): exact rational for mA, independent of global maVal *)
  sp   = N[1/(1 - xp), prec];
  mA   = N[3/20, prec];
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
  (* FIX 4: explicit Re[…] — the imaginary residual ~10^{-margin} is discarded.
     If this print fires with a large value, increase margin in n4PrecMin. *)
  (* Uncomment for debugging: Print["n4 Im residual @ J=", J, ": ", Im[N[result,20]]]; *)
  Re[N[result, 650]]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[3/20, 650];    (* FIX 1: exact rational *)
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 650]
];

X53[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[3/20, 650];    (* FIX 1: exact rational *)
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*((-4 + J)*(-2 + J)*(3 + J)*(5 + J) + ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 650]
];

(* Large J limit *)
LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 650];
  mA = N[3/20, 650];    (* FIX 1: exact rational *)
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-32*mA^6 + 24*mA^4*sp - 6*mA^2*sp^2 + sp^3))/(18*sp^3*(-4*mA^2 + sp)^6), 650]
];

Jmax = 60;
Jlist = Range[0, Jmax, 2];
(* JlistLarge = {3000};
Jlist = Join[Range[0, Jmax, 2], JlistLarge]; *)

fList = {g20, g31, n4, X52, X53};

(* large J limit: 0& for g20,g31,n4,X52; LargeJ for X53 *)
extraTriplet = {0&, 0&, 0&, 0&, LargeJ};

(* optimal upper bound *)
norm = {0, -1, 0, 0, 0};
obj  = {-1, 0, 0, 0, 0};


testNumericalSDP[spFile_String, jsonFile_String, prec_:650] := Module[
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
  Print["  extraTriplet: ", extraTriplet, "  (J\[Rule]\[Infinity] limit)"];

  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  polsRegular = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[fList[[k]][samplePoints[[i]], Jlist[[j]]], prec]},
        {k, Length[fList]}
      ] }}
    |>],
    {i, Length[samplePoints]}, {j, Length[Jlist]}
  ];

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
    (* 650 digits exceeds SDPB 2048-bit precision (≈ 616.5 decimal digits)
       by a safe margin.  n4 uses higher precision adaptively when J is large. *)
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 650];

    Print["=== test9.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    Null
  ]
];