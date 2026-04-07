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

f1[x_?NumericQ, J_?IntegerQ] := (1 + x)^2;
f2[x_?NumericQ, J_?IntegerQ] := ( 1 + x )*( 3 - 2*J*( J + 1 ));
f3[x_?NumericQ, J_?IntegerQ] := 2 * J * ( J + 1 )*( J*( J + 1 ) - 8 );

Jmax = 40;
Jlist = Range[0, Jmax, 2];  (* J = 0,2,...,40 *)

testNumericalSDP[jsonFile_, prec_:200] := Module[
  {samplePoints, sampleScalings, pols, norm, obj},

  samplePoints   = SetPrecision[{0.10340307765328547223630314082809872663649381367363570975866784459698950504846341436922826311233950415189685238409396227869590051487073360626174297880103181595478543065121171494081853010904892646771338, 0.97916391245035555397176210290443067298734173579221772588582301539489572733566938662472014152406657682245446430494243569601717924997604591473186319949654140608234966937487306130555629123626373267927949, 3.1450444531112909997984044773046581892845330789724131607992431051810113773940557833026483070649486883573151970746975847192422203679264870466595378300400920362599311630852176696089772742812014549423457}, prec];
  sampleScalings = SetPrecision[{0.90176341953323843798759045259032131659279233414342432294447471521607637811482288681937794999943292682649086768222789014656104058902714922501490993519844845467193976549344637063847467307127543462916940, 0.37562502300414320000623464943592436927747038445337677317825463674907535940473614239996104636093545947398985152128846241183359713001975271219984104672033251785836160106566703709943481372689361805633954, 0.043065009630614498560144075965094246234665231509999232187616552190861984641783256918190761425941162230338792571115907030330011934483910316041580807782687043748534787062965609301934098806334671855477152}, prec];

  (* One block per sample point. Each encodes a degree-0 (constant) polynomial
     constraint  a·f1(xᵢ) + b·f2(xᵢ) ≥ 0.

     Polynomials nesting:  {{{ {f1(xᵢ)}, {f2(xᵢ)} }}}
       {  }  ← Mathematica level 1 → JSON level 1 (column list, m_j=1)
        {  } ← Mathematica level 2 → JSON level 2 (row within column, m_j=1)
         {  }← Mathematica level 3 → JSON level 3 (polynomial vector, 2 elements)
       {f1(xᵢ)} ← Mathematica level 4 → JSON level 4 (degree-0 coefficient list)
       {f2(xᵢ)} ← Mathematica level 4 → JSON level 4 (degree-0 coefficient list)

     Verified against pmp.json which also has depth 4 for its polynomial strings. *)
  pols = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      "polynomials"    -> {{{               (* THREE wrapping pairs — NOT four *)
        {f1[ samplePoints[[i]], Jlist[[j]] ]},            (* Q^0: degree-0 coeff list, len=1 *)
        {f2[ samplePoints[[i]], Jlist[[j]] ]},             (* Q^1: degree-0 coeff list, len=1 *)
        {f3[ samplePoints[[i]], Jlist[[j]] ]}              (* Q^2: degree-0 coeff list, len=1 *)
      }}}
    |>],
    {i, Length[samplePoints]}, {j, Length[Jlist]}
  ];

  norm = {0, 1, 0};
  obj  = {-1, 0, 0};

  WritePmpJsonNumerical[jsonFile, SDP[obj, norm, Flatten[pols]], prec]
];

testNumericalSDP["numeric_pmp.json", 200];