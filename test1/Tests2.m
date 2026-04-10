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
   Define the n constraint functions as a list fVec = {f1, f2, …, fn}.
   Each function must accept a single numeric argument and return a
   numeric value.  No symbolic or polynomial form is required — only
   numeric evaluation at the sample points is ever used.

   Also set obj (objective vector b) and norm (normalisation vector n),
   both of length n = Length[fVec], such that:
     Maximise   b · y
     subject to n · y = 1   and   F(x) = sum_k fk(x)*yk >= 0 for all x

   EXAMPLE  (n = 2):
     Maximise -y₂  subject to  (1+x^4)·y₁ + (x^4/12+x^2)·y₂ ≥ 0
     obj  = {0, -1}  →  b·y = -y₂
     norm = {1,  0}  →  n·y = y₁ = 1
   ================================================================ *)

f1[x_?NumericQ] := 1 + x^4;
f2[x_?NumericQ] := x^4/12 + x^2;
(* Add further definitions here for n > 2, e.g.:
   f3[x_?NumericQ] := …; *)

fVec = {f1, f2};   (* ← replace {f1, f2, …, fn} to match your problem *)
norm = {1, 0};     (* ← n-vector: normalisation constraint  n·y = 1   *)
obj  = {0, -1};    (* ← n-vector: objective vector b  (maximise b·y)  *)

(* ================================================================
   END OF PROBLEM-SPECIFIC SECTION
   ================================================================ *)


(* ----------------------------------------------------------------
   testNumericalSDP
   ----------------------------------------------------------------
   Read sample points from spFile (one number per line, blank lines
   and lines starting with "#" are ignored), build one constant-
   function SDP block per point, and write the PMP JSON to jsonFile.

   The positivity constraint in each block is:
     fVec[[1]](xᵢ)·y₁ + fVec[[2]](xᵢ)·y₂ + … + fVec[[n]](xᵢ)·yₙ ≥ 0

   Only numeric evaluations of fVec[[k]] at xᵢ are used; the
   functions need not have any polynomial or closed-form structure.

   sampleScalings are derived as Exp[-xᵢ], the value of the
   DampedRational[1,{},1/E,x] prefactor at each xᵢ.  This formula
   is independent of the dimension n of y.
   ---------------------------------------------------------------- *)
testNumericalSDP[spFile_String, jsonFile_String, prec_:200] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, pols},

  (* --- Read and parse sampling_points.txt --- *)
  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " sample points from ", spFile];
  Print["  n = ", Length[fVec], " functional component(s)"];
  Print["  Points: ", samplePoints];

  (* --- Derive scalings from DampedRational[1,{},1/E,x] = e^{-x}.     --- *)
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  (* --- Build one constant-function SDP block per sample point.

     Each block encodes the constraint F(xᵢ) ≥ 0, i.e.
       fVec[[1]](xᵢ)·y₁ + … + fVec[[n]](xᵢ)·yₙ ≥ 0.

     Polynomials nesting:  {{ Table[{fVec[[k]][xᵢ]}, {k,n}] }}
       {{ … }}  ← levels 1 and 2 (column list / row within column, each size 1)
       Table[…] ← level 3: polynomial vector; n entries, one per component k
       {fVec[[k]][xᵢ]} ← level 4: degree-0 coefficient list (exactly 1 element)

     For n=2 this produces {{{ {f1(xᵢ)}, {f2(xᵢ)} }}} — identical to the
     original hardcoded form, verified against pmp.json. --- *)
  pols = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      (* TWO outer braces = JSON levels 1 & 2.
         Table gives level-3 list: { {f1(xᵢ)}, …, {fn(xᵢ)} }.
         Each inner {…} is the level-4 coefficient list (1 element, degree-0). *)
      "polynomials"    -> {{ Table[
        {SetPrecision[fVec[[k]][samplePoints[[i]]], prec]},
        {k, Length[fVec]}
      ] }}
    |>],
    {i, Length[samplePoints]}
  ];

  WritePmpJsonNumerical[jsonFile, SDP[obj, norm, pols], prec];
  Print["Wrote PMP JSON to ", jsonFile]
];

(* ----------------------------------------------------------------
   Command-line entry point
   ----------------------------------------------------------------
   USAGE:
     wolframscript -file Tests2.m <sample_points.txt> [output.json] [prec]

   ARGUMENTS:
     sample_points.txt   required  path to file with one x per line
     output.json         optional  output path (default: numeric_pmp.json)
     prec                optional  decimal digit precision (default: 200)

   EXAMPLE:
     wolframscript -file Tests2.m sampling_points.txt numeric_pmp.json 200

   NOTE ON ARGUMENT PARSING
     $ScriptCommandLine = {scriptname, arg1, arg2, ...}
     It is set identically by both  wolframscript -file  and  math -script,
     so no flag-hunting in $CommandLine is needed.
     When the file is loaded with  <<  as a library, $ScriptCommandLine
     is either empty or contains only the parent script's name, so
     Length[$ScriptCommandLine] < 2 and the block is skipped cleanly.
   ---------------------------------------------------------------- *)
Module[{myArgs, spFile, jsonFile, prec},

  (* Rest drops the leading script filename, leaving only user args. *)
  myArgs = If[Length[$ScriptCommandLine] >= 2,
    Rest[$ScriptCommandLine],
    {}
  ];

  If[Length[myArgs] >= 1,

    spFile   = myArgs[[1]];
    jsonFile = If[Length[myArgs] >= 2, myArgs[[2]], "numeric_pmp.json"];
    prec     = If[Length[myArgs] >= 3, ToExpression[myArgs[[3]]], 200];

    Print["=== Tests2.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    (* Loaded as library (<<) or called with no arguments — do nothing. *)
    Null
  ]
];