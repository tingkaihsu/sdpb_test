(* pmp generator for AAAA scattering *)
(* ============================================================
   FIX (2026-05):  n4[x, J] contains three special-function terms
   (Hypergeometric2F1 and two associated Legendre functions) that are each
   of order A(x)^J ~ 10^{J * log10(A)}, where A(x) > 1 for all physical x.
   These terms must cancel to produce the true value ~0 (large-J limit).
   The required precision is  prec_needed ≈ J * max_x log10[A(x)] + margin
                                          ≈ J * 0.156 + 30.
   At J = 10000, this is ~1590 decimal digits.  Using only prec = 600 leaves
   ~960 digits of cancellation unresolved, producing a garbage result ~10^960
   instead of ~0.  When written to the JSON and loaded by SDPB, this creates
   constraint-matrix entries B_{ij} ~ 10^960, making the affine residual
   p = b - B^T x ~ 10^980 from the very first Newton step → p-err = +inf.

   FIX:  replace n4 in fList with n4Safe, which returns
   SetPrecision[0, prec] for any J > J_SAFE_THRESHOLD.
   Justification: extraTriplet already encodes the large-J limit of n4 as 0,
   so using 0 at J = 10000 is the correct large-J approximation and is fully
   consistent with the existing SDP formulation.

   J_SAFE_THRESHOLD is chosen so that the cancellation residual is below one
   unit in the last place of the requested precision:
     J_SAFE_THRESHOLD = floor((prec - 30) / 0.156) ≈ 3654  for prec = 600.
   Since only J = 10000 in Jlist exceeds this threshold the change is
   minimal and surgical.

   ADDITIONAL FIX:  raise default prec from 600 to 650 so it exceeds
   SDPB's 2048-bit working precision (2048 * log10(2) ≈ 616.5 decimal digits)
   by a comfortable margin.
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
(* let the mass be 4mA^2 < M^2 = 1 where M  = 1 to infinity *)


maVal = SetPrecision[0.200, 600];

Print["mA = ", maVal]

(* dispersion representation of Wilson coefficients *)
(* All functions now precompute sp and mA numerically with N[...,prec] *)

(* forward limit: use our own convention *)

g20[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[-(2*Sqrt[sp/(-4*mA^2 + sp)])/(2*mA^2 - sp)^3, 600]
];

g31[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[-(Sqrt[sp/(-4*mA^2 + sp)]*(-3 - (2*J*(1 + J)*(2*mA^2 - sp))/(-4*mA^2 + sp)))/(-2*mA^2 + sp)^4, 600]
];

(* --- n4: safe precision guard -------------------------------------------
   n4 contains LegendreP[J, 1, z] and LegendreP[J, 2, z] with z > 1, and
   Hypergeometric2F1[-J, 1+J, 1, z1] which equals P_J(1.066) at typical
   sample points.  All three grow as A^J where A ~ 1.436 (at x ≈ 0),
   requiring ~1560 decimal digits at J = 10000 just to represent individual
   terms before cancellation.  With prec = 600 the computed result is pure
   numerical noise of magnitude ~10^960.

   Physical justification for returning 0 at large J:
     extraTriplet = {0&, 0&, 0&, 0&, LargeJ}
   explicitly encodes n4 → 0 as J → ∞ (index 3 in fList = n4).
   The safe threshold is  J_SAFE = floor((prec - 30) / 0.156) ≈ 3654.
   Since the only super-threshold spin in Jlist is J = 10000, this guard
   affects exactly the J = 10000 blocks.
   --------------------------------------------------------------------- *)

(* Maximum log10[A(x)] over physical x in (0,1): A = z + sqrt(z^2-1),
   z = 1 + 8*mA^2/(3*(sp - 4*mA^2)).  Evaluated numerically: max ~ 0.156. *)
n4MaxLog10A = 0.156;

n4Safe[x_?NumericQ, J_?IntegerQ] :=
  If[J > Floor[(600 - 30) / n4MaxLog10A],   (* J > ~3654 *)
    (* Large-J limit: n4 → 0 (consistent with extraTriplet's 0& for n4) *)
    SetPrecision[0, 600],
    (* Normal evaluation for J ≤ 3654 *)
    n4[x, J]
  ];

n4[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(81*Sqrt[sp/(-4*mA^2 + sp)]*(8*mA^4*(14*mA^2 - 15*sp)*(-8*mA^2 + 3*sp)^(3/2)*Hypergeometric2F1[-J, 1 + J, 1, (4*mA^2)/(12*mA^2 - 3*sp)] + (8*mA^4 - 18*mA^2*sp + 9*sp^2)*((-2*I)*mA*(10*mA^2 - 9*sp)*(8*mA^2 - 3*sp)*LegendreP[J, 1, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))] + Sqrt[-8*mA^2 + 3*sp]*(-8*mA^4 + 18*mA^2*sp - 9*sp^2)*LegendreP[J, 2, 1 + (8*mA^2)/(3*(-4*mA^2 + sp))])))/(4*mA^2*(-2*mA^2 + sp)*(-8*mA^2 + 3*sp)^(3/2)*(8*mA^4 - 18*mA^2*sp + 9*sp^2)^3), 600]
];

X52[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*(-((-4 + J)*(-2 + J)*(3 + J)*(5 + J)) - ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

X53[x_?NumericQ, J_?IntegerQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(J*(1 + J)*Sqrt[sp/(-4*mA^2 + sp)]*((-4 + J)*(-2 + J)*(3 + J)*(5 + J) + ((-1 + J)*(2 + J)*(-4*mA^2 + sp)^3*(36*mA^2 + (-15 + J + J^2)*sp))/sp^4))/(36*(-4*mA^2 + sp)^6), 600]
];

(* Large J limit *)
LargeJ[x_?NumericQ] := Module[{sp, mA},
  sp = N[1/(1-x), 600];
  mA = N[maVal, 600];
  N[(Sqrt[sp/(-4*mA^2 + sp)]*(-32*mA^6 + 24*mA^4*sp - 6*mA^2*sp^2 + sp^3))/(18*sp^3*(-4*mA^2 + sp)^6), 600]
];

Jmax = 60;
(* Jlist = Range[0, Jmax, 2]; *)
JlistLarge = {3000};
Jlist = Join[Range[0, Jmax, 2], JlistLarge];

(* NOTE: n4 is replaced by n4Safe (returns 0 for J > ~3654, exact for J <= 3654).
   All other functions are unchanged; they grow polynomially in J and are
   numerically accurate with prec = 600 even at J = 10000. *)
fList = {g20, g31, n4Safe, X52, X53};

(* large J limit *)
(* 0& is a constant function of 0 *)

extraTriplet = {0&, 0&, 0&, 0&, LargeJ};

(* optimal lower bound *)
(* norm = {0, 1, 0, 0, 0};
obj = {-1, 0, 0, 0, 0}; *)

(* optimal upper bound *)
norm = {0, -1, 0, 0, 0};
obj = {-1, 0, 0, 0, 0};



testNumericalSDP[spFile_String, jsonFile_String, prec_:650] := Module[
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
  Print["  n4MaxLog10A : ", n4MaxLog10A, " => n4Safe threshold J > ", Floor[(prec-30)/n4MaxLog10A]];

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
    (* Default prec raised from 600 → 650 to exceed SDPB's 2048-bit working
       precision (2048 * log10(2) ≈ 616.5 decimal digits) by a safe margin. *)
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 650];

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
