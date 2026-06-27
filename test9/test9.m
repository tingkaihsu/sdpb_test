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
LaunchKernels[];

(* ================================================================
   PROBLEM-SPECIFIC SECTION  ← edit here for your problem
   ----------------------------------------------------------------
   Functions fk[x, J] depend on BOTH x (continuous, discretised)
   AND J (discrete spin, constrained EXACTLY for J in Jlist).
   Sampling is done over x; J constraints are imposed exactly.

   functionalVector[x,J] collects {g2, g3} and the selected null constraints into one vector.
   sp and mA are assigned once inside functionalVector for each (x,J) evaluation.

   Jmax and Jlist define the discrete spin sum (even spins only).

   largeJVector[x] is the J → ∞ limiting constraint vector:
     As J → ∞, divide fk[x,J] by the leading power of J (here J^4):
       all listed components vanish except the x82 leading term
       largeJVector[x] stores that x82 leading coefficient in the x82 slot
       largeJVector keeps this vector aligned with functionalVector
     Only the x82 component is nonzero in the large-J vector.
     This enforces the sampled large-J limiting inequality in the selected dual basis.
   ================================================================ *)

(* 4mA^2 < M^2 *)
mAval = N[1/5, 1000];
Print["mA = ", mAval];

(* Bootstrap target entries: {g2, g3}.  Null constraints are appended below. *)
cList[sp_, mA_, J_] := {
  -((2 Sqrt[sp/(-4 mA^2+sp)])/(2 mA^2-sp)^3),
  (Sqrt[sp/(-4 mA^2+sp)] (-3-(2 J (1+J) (2 mA^2-sp))/(-4 mA^2+sp)))/(-2 mA^2+sp)^4
};

(* Add/remove null constraints here, then update nulllist by label. *)
nulllist = {6, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];

Nlist[n_, sp_, mA_, J_] := {-((243 Sqrt[-(sp/(4 mA^2-sp))] (6 (4 mA^2-sp) Sqrt[(8 mA^4-3 mA^2 sp)/(-4 mA^2+sp)^2] LegendreP[J,1,1+(8 mA^2)/(3 (-4 mA^2+sp))]+(4 mA^2-3 sp) LegendreP[J,2,1+(8 mA^2)/(3 (-4 mA^2+sp))]))/(4 mA^2 (4 mA^2-3 sp)^4 (8 mA^2-3 sp))),((J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-4+J) (-2+J) (3+J) (5+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (36 mA^2+(-15+J+J^2) sp))/sp^4))/(36 (-4 mA^2+sp)^6)),(J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J))-((-1+J) (2+J) (-4 mA^2+sp)^3 (-2304 mA^4+1152 mA^2 sp+(-72+J (1+J) (-18+J+J^2)) sp^2))/sp^5))/(576 (-4 mA^2+sp)^7),(J (1+J) Sqrt[sp/(-4 mA^2+sp)] (-(((-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (-47+J+J^2))/(-4 mA^2+sp)^8)+((-1+J) (2+J) (((-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) sp^3)/(4 mA^2-sp)^5+3600/(-4 mA^2+sp)^2))/sp^6))/14400,((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-6+J) (-4+J) (5+J) (7+J)-((-1+J) (2+J) (-4 mA^2+sp)^4 (64 mA^2+(-28+J+J^2) sp))/sp^5))/(576 (-4 mA^2+sp)^8),(J (1+J) Sqrt[sp/(-4 mA^2+sp)] (((-8+J) (-6+J) (-4+J) (-2+J) (3+J) (5+J) (7+J) (9+J) (-38+J+J^2))/(4 mA^2-sp)^9+((-1+J) (2+J) (-(((-5+J) (-4+J) (-3+J) (-2+J) (3+J) (4+J) (5+J) (6+J) sp^4)/(-4 mA^2+sp)^6)+129600/(-4 mA^2+sp)^2))/sp^7))/518400,(((-2+J) J (1+J) (3+J) Sqrt[sp/(-4 mA^2+sp)] ((-8+J) (-6+J) (-4+J) (5+J) (7+J) (9+J)+((-2+J+J^2) (-4 mA^2+sp)^4 (6400 mA^4-3200 mA^2 sp+(160-J (1+J) (-32+J+J^2)) sp^2))/sp^6))/(14400 (-4 mA^2+sp)^9))}[[n+1]];
largeJlist[sp_, mA_] := {0, 0, 0, 0, 0, 0, 0, -((Sqrt[sp/(-4 mA^2+sp)] (-32 mA^6+24 mA^4 sp-6 mA^2 sp^2+sp^3))/(259200 sp^3 (-4 mA^2+sp)^9)), 0};

functionalVector[x_?NumericQ, J_?IntegerQ, prec_:1000] := Module[{sp, mA},
  sp = N[1/(1-x), prec];
  mA = N[mAval, prec];
  N[Join[cList[sp, mA, J], Table[Nlist[n,sp,mA,J],{n,0,nulllist[[1]]}] ], prec]
];


largeJVector[x_?NumericQ, prec_:1000] := Module[{sp, mA},
  sp = N[1/(1-x), prec];
  mA = N[mAval, prec];
  N[largeJlist[sp, mA], prec]
];

functionalDimension = 2 + Length[list0];

constantPolynomialVector[values_List, prec_] := {{SetPrecision[{#}, prec] & /@ values}};

samplePmp[polynomialValues_List, xi_, sampleScaling_, prec_] :=
  NumericalPositiveMatrixWithPrefactor[<|
    "prefactor"      -> DampedRational[1, {}, 1/E, x],
    "samplePoints"   -> {xi},
    "sampleScalings" -> {sampleScaling},
    "polynomials"    -> constantPolynomialVector[polynomialValues, prec]
  |>];

regularFunctionalValues[xi_, J_, prec_] := functionalVector[xi, J, prec];

largeJFunctionalValues[xi_, prec_] := largeJVector[xi, prec];

regularPmp[xi_, sampleScaling_, J_, prec_] :=
  samplePmp[regularFunctionalValues[xi, J, prec], xi, sampleScaling, prec];

largeJPmp[xi_, sampleScaling_, prec_] :=
  samplePmp[largeJFunctionalValues[xi, prec], xi, sampleScaling, prec];

buildNumericalPMPs[samplePoints_List, sampleScalings_List, prec_] := Module[
  {regularPMPs, largeJPMPs},

  regularPMPs = Flatten[{
    Flatten[ParallelTable[regularPmp[samplePoints[[i]], sampleScalings[[i]], J, prec],{i, Length[samplePoints]}, {J, 0, 1000, 2}]],
    Flatten[ParallelTable[regularPmp[samplePoints[[i]], sampleScalings[[i]], J, prec],{i, Length[samplePoints]}, {J, 1500, 5000, 100}]],
    Flatten[ParallelTable[regularPmp[samplePoints[[i]], sampleScalings[[i]], J, prec],{i, Length[samplePoints]}, {J, 6000, 20000, 500}]],
    Flatten[ParallelTable[regularPmp[samplePoints[[i]], sampleScalings[[i]], J, prec],{i, Length[samplePoints]}, {J, 20000, 50000, 2000}]]
  },1];

  largeJPMPs = Table[
    largeJPmp[samplePoints[[i]], sampleScalings[[i]], prec],
    {i, Length[samplePoints]}
  ];

  Join[Flatten[regularPMPs], largeJPMPs]
];

(* test9 targets g3/g2: objective on g2, normalization on g3. *)
norm = -1 * Flatten[{{0,1},list0}];
obj  = -1 * Flatten[{{1,0},list0}];

Print["size of nomr = ", Length[norm]];
Print["size of obj = ", Length[obj]];

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

   For each xi an EXTRA block is created from largeJVector,
   enforcing the J→∞ limit constraint:
     largeJVector[xi][[1]]·y1 + ... + largeJVector[xi][[K]]·yK ≥ 0

   Total blocks = Length[samplePoints] × (Length[Jlist] + 1).

   sampleScalings are Exp[-xi], the value of the prefactor
   DampedRational[1,{},1/E,x] at each sample point xi.
   ---------------------------------------------------------------- *)
ToPMP[spFile_String, jsonFile_String, prec_:1000] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, pols},

  (* --- Read and parse sampling_points.txt --- *)
  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " x sample points from ", spFile];
  
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];
  pols = buildNumericalPMPs[samplePoints, sampleScalings, prec];

  WritePmpJsonNumerical[jsonFile, SDP[obj, norm, pols], prec];
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

    Print["=== test9.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    ToPMP[spFile, jsonFile, prec];
    Quit[0],

    Null
  ]
];