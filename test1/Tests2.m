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

f1[x_?NumericQ] := 1 + x^4;
f2[x_?NumericQ] := x^4/12 + x^2;

(* ----------------------------------------------------------------
   testNumericalSDP
   ----------------------------------------------------------------
   Read sample points from spFile (one number per line, blank lines
   and lines starting with "#" are ignored), build one constant-
   function SDP block per point, and write the PMP JSON to jsonFile.

   sampleScalings are derived as Exp[-xᵢ], i.e. the value of the
   DampedRational[1,{},1/E,x] prefactor at each sample point.
   This is the correct systematic choice: the scaling for block i
   should equal the prefactor evaluated at xᵢ.
   ---------------------------------------------------------------- *)
testNumericalSDP[spFile_String, jsonFile_String, prec_:200] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, pols, norm, obj},

  (* --- Read and parse sample_points.txt --- *)
  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " sample points from ", spFile];
  Print["  Points: ", samplePoints];

  (* --- Derive scalings from the prefactor DampedRational[1,{},1/E,x]
         evaluated at each xᵢ: value = (1/E)^xᵢ = Exp[-xᵢ].           --- *)
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  (* --- Build one constant-function block per sample point.

     One block per sample point. Each encodes a degree-0 (constant) polynomial
     constraint  a·f1(xᵢ) + b·f2(xᵢ) ≥ 0.

     Polynomials nesting:  {{{ {f1(xᵢ)}, {f2(xᵢ)} }}}
       {  }  ← Mathematica level 1 → JSON level 1 (column list, m_j=1)
        {  } ← Mathematica level 2 → JSON level 2 (row within column, m_j=1)
         {  }← Mathematica level 3 → JSON level 3 (polynomial vector, 2 elements)
       {f1(xᵢ)} ← Mathematica level 4 → JSON level 4 (degree-0 coefficient list)
       {f2(xᵢ)} ← Mathematica level 4 → JSON level 4 (degree-0 coefficient list)

     Verified against pmp.json which also has depth 4 for its polynomial strings. --- *)
  pols = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials"    -> {{{               (* THREE wrapping pairs — NOT four *)
        {f1[samplePoints[[i]]]},            (* Q^0: degree-0 coeff list, len=1 *)
        {f2[samplePoints[[i]]]}             (* Q^1: degree-0 coeff list, len=1 *)
      }}}
    |>],
    {i, Length[samplePoints]}
  ];

  norm = {1, 0};
  obj  = {0, -1};

  WritePmpJsonNumerical[jsonFile, SDP[obj, norm, pols], prec];
  Print["Wrote PMP JSON to ", jsonFile]
];

(* ----------------------------------------------------------------
   Command-line entry point
   ----------------------------------------------------------------
   USAGE:
     wolframscript -file Tests2.m <sample_points.txt> [output.json] [prec]
     math          -script Tests2.m <sample_points.txt> [output.json] [prec]

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