(* New head for purely numerical data *)
ClearAll[NumericalPositiveMatrixWithPrefactor];

(* Helper: recursively stringify all numeric leaves *)
toJsonNestedNumberArray[expr_, prec_] := expr /. n_?NumericQ :> toJsonNumber[n, prec];

(*
  Design: the "constant-function" approach.

  For each sample point xᵢ we create ONE PositiveMatrixWithPrefactor block whose
  polynomial matrix entries are CONSTANT (degree-0) polynomials:
      W⁰ⱼ(x) = f1(xᵢ),   W¹ⱼ(x) = f2(xᵢ)

  This gives n independent 1×1 PSD constraints solved simultaneously:
      a·f1(xᵢ) + b·f2(xᵢ) ≥ 0   for each i = 1..n

  SDPB manual (Section 3.1, p.4) JSON nesting for "polynomials":
    Level 1: [ ]   — column list          (m_j columns; here m_j=1)
    Level 2:  [ ]  — row within column    (m_j rows;    here m_j=1)
    Level 3:   [ ] — polynomial vector    [Q^0, Q^1]  (N+1=2 entries for N=1)
    Level 4:    [ "c0", ..., "cd" ]       coefficient list of each polynomial

  For degree-0 (constant) polynomials dⱼ=0, so each coefficient list has 1 element:
    Level 4: [ "f1(xᵢ)" ]   and   [ "f2(xᵢ)" ]

  CORRECT JSON output for block i:
    "polynomials": [              ← level 1 (len=1)
        [                         ← level 2 (len=1)
            [                     ← level 3 (len=2): polynomial vector
                ["f1(xᵢ)"],       ← level 4 (len=1): Q^0 coefficient list
                ["f2(xᵢ)"]        ← level 4 (len=1): Q^1 coefficient list
            ]
        ]
    ]

  In Mathematica this requires exactly THREE wrapping brace pairs before the
  coefficient lists:
      {{{  {f1(xᵢ)}, {f2(xᵢ)}  }}}
       ↑1   ↑2        ↑↑ depth-3 elements = degree-0 coefficient lists

  FOUR wrapping pairs (the previous bug):
      {{{{  {f1(xᵢ)}, {f2(xᵢ)}  }}}}
  produced depth 5, pushing coefficient lists one level too deep. The parser,
  while in the Json_Float_Parser state (expecting float strings at depth 4),
  received another '[' (array start), triggering:
      "Not implemented: function 'bool json_start_array()' in class:
       '17Json_Float_ParserIN2El8BigFloatEE'"
*)

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

(* ================================================================
   PROBLEM-SPECIFIC SECTION  ← edit here for your problem
   ----------------------------------------------------------------
   Functions fk[x, J] depend on BOTH x (continuous, discretised)
   AND J (discrete spin, constrained EXACTLY for J in Jlist).
   Sampling is done over x; J constraints are imposed exactly.

   fList collects the functions for the generalised Table loop.
   Each fList[[k]][x, J] must return a numeric value.

   Jmax and Jlist define the discrete spin sum (even spins only).

   extraTriplet is the J → ∞ limiting constraint vector:
     As J → ∞, divide fk[x,J] by the leading power of J (here J^4):
       f1[x,J] / J^4 → 0        (f1 is J-independent)
       f2[x,J] / J^4 → 0        (f2 ~ J^2, subdominant)
       f3[x,J] / J^4 → 2        (f3 ~ 2·J^4, leading term)
     So extraTriplet = {0, 0, 2}.
     This enforces  0·y1 + 0·y2 + 2·y3 ≥ 0  in the J→∞ limit.
   ================================================================ *)
m1 = N[1/2, 600];
J1 = 0;
J2 = 2;
mgap = N[166/100, 600];

(* heavy-sum *)
g0[x_?NumericQ, J_?IntegerQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[2/s, 600]
];

x10[x_?NumericQ, J_?IntegerQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[(2-J (7+J))/(2 s^2), 600]
];

x20[x_?NumericQ, J_?IntegerQ] :=  Module[{s},
  s = N[mgap/(1-x), 600];
  N[(20-J (7+J) (-13+J (7+J)))/(20 s^3), 600]
];

x30[x_?NumericQ, J_?IntegerQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 s^4)), 600]
];

x40[x_?NumericQ, J_?IntegerQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[1/s^5-1/(10080 s^5)(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))), 600]
];

x41[x_?NumericQ, J_?IntegerQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[-((J (7+J) (-23+J (7+J)))/(20 s^5)), 600]
];

x50[x_?NumericQ, J_?IntegerQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[1/s^6-1/(403200 s^6)(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)), 600]
];

largeJ[x_?NumericQ] := Module[{s},
  s = N[mgap/(1-x), 600];
  N[-(1/(403200 s^6)), 600]
];

(* first state resonance *)

g0fst[x_?NumericQ ] := Module[{s},
  s = N[m1, 600];
  N[2/s, 600]
];

x10fst[x_?NumericQ ] := Module[{s, J},
  s = N[m1, 600];
  J = J1;
  N[(2-J (7+J))/(2 s^2), 600]
];

x20fst[x_?NumericQ ] := Module[{s, J},
  s = N[m1, 600];
  J = J1;
  N[(20-J (7+J) (-13+J (7+J)))/(20 s^3), 600]
];

x30fst[x_?NumericQ ] := Module[{s, J},
  s = N[m1, 600];
  J = J1;
  N[-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 s^4)), 600]
];

x40fst[x_?NumericQ ] := Module[{s, J},
  s = N[m1, 600];
  J = J1;
  N[1/s^5-1/(10080 s^5)(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))), 600]
];

x41fst[x_?NumericQ ] := Module[{s, J},
  s = N[m1, 600];
  J = J1;
  N[-((J (7+J) (-23+J (7+J)))/(20 s^5)), 600]
];

x50fst[x_?NumericQ ] := Module[{s, J},
  s = N[m1, 600];
  J = J1;
  N[1/s^6-1/(403200 s^6)(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)), 600]
];

(* norm vector (snd state resonance) *)

g0snd = Module[{s},
  s = N[1, 600];
  N[2/s, 600]
];

x10snd = Module[{s, J},
  s = N[1, 600];
  J = J2;
  N[(2-J (7+J))/(2 s^2), 600]
];

x20snd = Module[{s, J},
  s = N[1, 600];
  J = J2;
  N[(20-J (7+J) (-13+J (7+J)))/(20 s^3), 600]
];

x30snd = Module[{s, J},
  s = N[1, 600];
  J = J2;
  N[-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 s^4)), 600]
];

x40snd = Module[{s, J},
  s = N[1, 600];
  J = J2;
  N[1/s^5-1/(10080 s^5)(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))), 600]
];

x41snd = Module[{s, J},
  s = N[1, 600];
  J = J2;
  N[-((J (7+J) (-23+J (7+J)))/(20 s^5)), 600]
];

x50snd = Module[{s, J},
  s = N[1, 600];
  J = J2;
  N[1/s^6-1/(403200 s^6)(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)), 600]
];

(* fList = {g0, n2, n4, x52, x62, x72};   (* one entry per component of y; must match Length[norm] *)
jList = {0&, 0&, 0&, 0&, 0&, largeJ};
ResList = {g0shft, n2shft, n4shft, x52shft, x62shft, x72shft}; *)

fList = {g0, x10, x20, x30, x40, x41, x50};
jList = {0&, 0&, 0&, 0&, 0&, 0&, largeJ};
ResList = {g0fst, x10fst, x20fst, x30fst, x40fst, x41fst, x50fst};

(* norm = {G0, N2, N4, X52, X62, X72}; *)
norm = {g0snd, x10snd, x20snd, x30snd, x40snd, x41snd, x50snd};

Jmax  = 60;
Jlist = Range[0, Jmax, 2];   (* J = 0, 2, 4, …, 40 — exact discrete constraints *)

(* obj  = {-1, 0, 0, 0, 0, 0};   objective: maximise -y1 = minimise y1 *)

obj = {-1, 0, 0, 0, 0, 0, 0};

(* ================================================================
   END OF PROBLEM-SPECIFIC SECTION
   ================================================================ *)


(* ----------------------------------------------------------------
   testNumericalSDP
   ----------------------------------------------------------------
   Read x sample points from spFile (one number per line; blank
   lines and "#"-comment lines are ignored).

   For each (xi, Jj) pair a REGULAR block is created enforcing:
     f1(xi,Jj)·y1 + f2(xi,Jj)·y2 + f3(xi,Jj)·y3 ≥ 0

   For each xi an EXTRA block is created from extraTriplet,
   enforcing the J→∞ limit constraint:
     extraTriplet[[1]]·y1 + extraTriplet[[2]]·y2 + extraTriplet[[3]]·y3 ≥ 0

   Total blocks = Length[samplePoints] × (Length[Jlist] + 1).

   sampleScalings are Exp[-xi], the value of the prefactor
   DampedRational[1,{},1/E,x] at each sample point xi.
   ---------------------------------------------------------------- *)
testNumericalSDP[spFile_String, jsonFile_String, prec_:1000] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, polsRegular, polsExtra},

  (* --- Read and parse sampling_points.txt --- *)
  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " x sample points from ", spFile];
  (* Print["  x-points    : ", samplePoints];
  Print["  J-values    : ", Jlist, "  (", Length[Jlist], " spins, exact)"];
  Print["  extraTriplet: ", extraTriplet, "  (J\[Rule]\[Infinity] limit)"]; *)

  (* Scalings = prefactor DampedRational[1,{},1/E,x] evaluated at xi = e^{-xi} *)
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  (* --- Regular blocks: one per (xi, Jj) pair.
     Polynomials nesting:  {{ Table[{fk(xi,Jj)}, {k,3}] }}
       {{ ... }}  ← JSON levels 1 and 2  (column list / row, each size 1)
       Table[...] ← JSON level 3: polynomial vector, one entry per fList[[k]]
       {fk(xi,Jj)} ← JSON level 4: degree-0 coefficient list (1 element) --- *)

  pols = Table[
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

  (* --- Extra blocks: one per xi, encoding the J→∞ limit constraint.
     extraTriplet = {0, 0, 2} is x-independent (J^4 leading coefficient),
     so every extra block carries the same polynomial values.
     Each block is still associated with a distinct xi for SDPB bookkeeping. --- *)
  polsLargeJ = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[jList[[k]][samplePoints[[i]]], prec]},
        {k, Length[jList]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  pols1Res = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[ResList[[k]][samplePoints[[i]]], prec]},
        {k, Length[ResList]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  (* Flatten polsRegular (2D Table → flat list) and append polsExtra (already flat) *)
  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Join[Flatten[pols], polsLargeJ, pols1Res]],
    prec
  ];
  Print["Wrote PMP JSON to ", jsonFile]
];


(* ----------------------------------------------------------------
   Command-line entry point
   ----------------------------------------------------------------
   USAGE:
     wolframscript -file g3_ExtremalEFT_2.m <sp_file> [output.json] [prec]

   ARGUMENTS:
     sp_file       required  path to sampling_points.txt (one x per line)
     output.json   optional  output path (default: numeric_pmp.json)
     prec          optional  decimal digit precision (default: 200)

   NOTE: $ScriptCommandLine = {scriptname, arg1, arg2, …} is populated
   identically by both wolframscript -file and math -script. Rest[] drops
   the script name leaving only user arguments. When loaded with << as a
   library, $ScriptCommandLine has length ≤ 1 and the block is skipped.
   ---------------------------------------------------------------- *)
Module[{myArgs, spFile, jsonFile, prec},

  myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

  If[Length[myArgs] >= 1,
    spFile   = myArgs[[1]];
    jsonFile = If[Length[myArgs] >= 2, myArgs[[2]], "numeric_pmp.json"];
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 1000];

    Print["=== g3_ExtremalEFT_2.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    Null
  ]
];