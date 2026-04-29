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

f1[x_?NumericQ, J_?IntegerQ] := 2*((1-x))^3;
f2[x_?NumericQ, J_?IntegerQ] := -(3 - 2*J*(J + 1))*((1-x))^4;


X52[x_?NumericQ, J_?IntegerQ] := -1/36*(J*(1 + J)*(150 + J*(1 + J)*(-43 + 2*J*(1 + J))));
X62[x_?NumericQ, J_?IntegerQ] := -1/288*((-3 + J)*J*(1 + J)*(4 + J)*(204 + J*(1 + J)*(-32 + J + J^2)));
X72[x_?NumericQ, J_?IntegerQ] := -1/14400*(J*(1 + J)*(246960 + J*(1 + J)*(-67908 + J*(1 + J)*(4916 + J*(1 + J)*(-155 + 2*J*(1 + J))))));
X82[x_?NumericQ, J_?IntegerQ] := -1/259200*(J*(1 + J)*(-6808320 + J*(1 + J)*(1906416 + J*(1 + J)*(-170976 + J*(1 + J)*(6568 + J*(1 + J)*(-124 + J + J^2))))));
X92[x_?NumericQ, J_?IntegerQ] := -1/25401600*(J*(1 + J)*(1015701120 + J*(1 + J)*(-306848736 + J*(1 + J)*(28977336 + J*(1 + J)*(-1293996 + J*(1 + J)*(30170 + J*(1 + J)*(-371 + 2*J*(1 + J))))))));
X102[x_?NumericQ, J_?IntegerQ] := -1/812851200*(J*(1 + J)*(-44242329600 + J*(1 + J)*(13817329920 + J*(1 + J)*(-1475388288 + J*(1 + J)*(74195472 + J*(1 + J)*(-2018816 + J*(1 + J)*(31080 + J*(1 + J)*(-264 + J + J^2))))))));


fList = {f1, f2, X52, X62, X72, X82, X92, X102};   (* one entry per component of y; must match Length[norm] *)

Jmax  = 40;
Jlist = Range[0, Jmax, 2];   (* J = 0, 2, 4, …, 40 — exact discrete constraints *)

(* Extremal-J limit: leading J^4 coefficient vector of {f1, f2, f3}.
   x-independent, so the same vector is used for all sample points xi. *)

extraTriplet = {0, 0, 0, 0, 0, 0, 0, -1/812851200};

norm = {0, 1, 0, 0, 0, 0, 0, 0};    (* normalisation: y2 = 1                *)
obj  = {-1, 0, 0, 0, 0, 0, 0, 0};   (* objective: maximise -y1 = minimise y1 *)

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
testNumericalSDP[spFile_String, jsonFile_String, prec_:200] := Module[
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

  (* --- Extra blocks: one per xi, encoding the J→∞ limit constraint.
     extraTriplet = {0, 0, 2} is x-independent (J^4 leading coefficient),
     so every extra block carries the same polynomial values.
     Each block is still associated with a distinct xi for SDPB bookkeeping. --- *)
  polsExtra = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials" -> {{ Table[
        {SetPrecision[extraTriplet[[k]], prec]},
        {k, Length[extraTriplet]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  Print["  Regular blocks : ", Length[samplePoints] * Length[Jlist]];
  Print["  Extra blocks   : ", Length[polsExtra], "  (J\[Rule]\[Infinity] constraint)"];
  Print["  Total blocks   : ", Length[samplePoints] * Length[Jlist] + Length[polsExtra]];

  (* Flatten polsRegular (2D Table → flat list) and append polsExtra (already flat) *)
  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Join[Flatten[polsRegular], polsExtra]],
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
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 200];

    Print["=== g3_ExtremalEFT_2.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    (* Loaded with << as a library — do nothing. *)
    Null
  ]
];