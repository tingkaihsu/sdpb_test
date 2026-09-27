(* ::Package:: *)

(* Exported from AAAA_sum_rules.original.nb. *)
(* All original Input cells are included in their original order. *)
(* Only the Linear Dependence block has been replaced by the optimized checker. *)

(* ===== Input cell 1 ===== *)
PartialWaveD[d_,J_,x_]:=Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-x)/2];


PartialWaveD[4,J,x]-LegendreP[J,x]//FullSimplify

(* ===== Input cell 2 ===== *)
 ((x-4*mA^2)^((10-3)/(2)))/(Sqrt[x])/.{mA->0}
 ((x-4*mA^2)^((10-3)/(2)))/(Sqrt[x])//FullSimplify

(* ===== Input cell 3 ===== *)
(* null constraint generator for AAAA scattering *)
(* ---------- coefficient helper ---------- *)

del[a_, b_] := delCoeff @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* AAAA scattering change the ansatz to be s-u symmetric *)
Mfwdlow[s_, t_, mA_, Nmax_Integer] :=
    Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s-2*mA^2)^ab[[1]]
            * (u-2*mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] /. {u -> 4*mA^2 - s - t};

(* ===== Input cell 4 ===== *)
Ker[s_,t_,s1_,s2_,k_, mA_]:=(Mfwdlow[s,t,mA,10])/((s-s1)*((s-s1)*(s-s2))^(k/2));


s1 = 2*mA^2;


s2 = 2*mA^2-t;

(* ===== Input cell 5 ===== *)
(* g0 + ... *)
-Residue[Ker[s,t,s1,s2,0,mA],{s,Infinity}]//FullSimplify

(* ===== Input cell 6 ===== *)
SumRule[s_,t_,s1_,s2_,k_,J_,mA_]:=((1)/((s-s1)*((s-s1)*(s-s2))^(k/2))-1/((4*mA^2-s-t-s1)*((4*mA^2-s-t-s1)*(4*mA^2-s-t-s2))^(k/2)))*(((-4 (mA)^(2)+s))^(7/2))/(Sqrt[s])*PartialWaveD[10,J,1+(2*t)/(s-4*mA^2)];

(* ===== Input cell 7 ===== *)
SeriesCoefficient[SumRule[x,t,s1,s2,0,J,mA],{t,0,0}]//FullSimplify

(* ===== Input cell 8 ===== *)
DblCtr[mA_, k_Integer, q_Integer, Nmax_Integer] := Residue[ Residue[ 1/(s*t) ( Mfwdlow[s, t, mA, Nmax]/(s^(k-q)*t^q) - Mfwdlow[t, s, mA, Nmax]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

(* ===== Input cell 9 ===== *)
(* Double Contour Integral *)
DblCtrSum[s_,k_,q_]:=Residue[((1)/(s*t)*((1)/(s^(k-q)*t^q)-(1)/(s^q*t^(k-q)))-(1)/((4*mA^2-s-t)*t)*((1)/((4*mA^2-s-t)^(k-q)*t^q)-(1)/((4*mA^2-s-t)^q*t^(k-q))))*(((-4 (mA)^(2)+s))^(7/2))/(Sqrt[s])*PartialWaveD[10,J,1+(2*t)/(s-4*mA^2)],{t,0}];

(* ===== Input cell 10 ===== *)
x10 = DblCtrSum[x,1,0]//FullSimplify;

(* ===== Input cell 11 ===== *)
x20 = DblCtrSum[x,2,0]//FullSimplify;

(* ===== Input cell 12 ===== *)
x30 = DblCtrSum[x,3,0]//FullSimplify;


x31 = DblCtrSum[x,3,1]//FullSimplify;

(* ===== Input cell 13 ===== *)
x40 = DblCtrSum[x,4,0]//FullSimplify;


x41 = DblCtrSum[x,4,1]//FullSimplify;

(* ===== Input cell 14 ===== *)
x50=DblCtrSum[x,5,0]//FullSimplify;


x51 = DblCtrSum[x,5,1];


x52 = DblCtrSum[x,5,2];

(* ===== Input cell 15 ===== *)
x60 = DblCtrSum[x,6,0];


x61 = DblCtrSum[x,6,1];


x62 = DblCtrSum[x,6,2];

(* ===== Input cell 16 ===== *)
x70 = DblCtrSum[x,7,0];


x71 = DblCtrSum[x,7,1];


x73 = DblCtrSum[x,7,3];

(* ===== Input cell 17 ===== *)
x80 = DblCtrSum[x,8,0];


x81 = DblCtrSum[x,8,1];


x83 = DblCtrSum[x,8,3];

(* ===== Input cell 18 ===== *)
x90 = DblCtrSum[x,9,0];


x91 = DblCtrSum[x,9,1];


x93 =DblCtrSum[x,9,3];
x94 = DblCtrSum[x,9,4];


x100 = DblCtrSum[x,10,0];
x101 = DblCtrSum[x,10,1];
x103 = DblCtrSum[x,10,3];
x104 = DblCtrSum[x,10,4];

x110 = DblCtrSum[x,11,0];
x111 = DblCtrSum[x,11,1];
x112 = DblCtrSum[x,11,2];
x113 = DblCtrSum[x,11,3];
x114 = DblCtrSum[x,11,4];
x115 = DblCtrSum[x,11,5];


(* ===== Input cell 19 ===== *)
ClearAll[largeJ];



largeJ[list_List,J_Symbol:J]:=Module[{p},p=Max[Exponent[Together@Simplify[#],J]&/@list];
Assuming[J>0,Limit[list/J^p,J->Infinity]]]

(* ===== Input cell 20 ===== *)
lst = {x10,x20,x30, x40, x41 ,x50, x51,x52,x60,x61,x62,x70,x71,x73,x80,x81,x83,x90,x91,x93,x94,x100,x101,x103,x104,x110,x111,x112,x113,x114,x115};



largeJ[lst]//FullSimplify



Length[lst]

(* ===== Input cell 21 ===== *)
lst//FullSimplify

(* ===== Input cell 22 ===== *)
(* Adaptive exact linear-dependence checker for null constraints. *)

ClearAll[
  HoldNames,
  DependenceAtoms,
  SafeRationalRules,
  EvaluationVector,
  ProbeMatrix,
  RankIncrementIndependentColumns,
  VerifyRelation,
  RelationToRule,
  CheckNullDependenceByProbing
];

SetAttributes[HoldNames, HoldAll];
HoldNames[syms_List] := HoldForm /@ Unevaluated[syms];

Options[DependenceAtoms] = {
  "ExtraAtoms" -> {},
  "AtomHeads" -> {delCoeff}
};

DependenceAtoms[exprs_, OptionsPattern[]] := Module[
  {heads, extra, held, found},
  heads = OptionValue["AtomHeads"];
  extra = OptionValue["ExtraAtoms"];
  held = Hold[exprs];
  found = DeleteDuplicates @ Cases[
      held,
      f_[___] /; MemberQ[heads, Unevaluated[f]],
      {0, Infinity}
    ];
  DeleteDuplicates @ Join[extra, found]
];

Options[SafeRationalRules] = {
  "IntegerRange" -> {-7, 7},
  "ExcludeZeroFor" -> {},
  "Seed" -> Automatic
};

SafeRationalRules[vars_List, OptionsPattern[]] := Module[
  {range, seed, excludeZero, values, pick},
  range = OptionValue["IntegerRange"];
  seed = OptionValue["Seed"];
  excludeZero = OptionValue["ExcludeZeroFor"];
  If[seed =!= Automatic, SeedRandom[seed]];
  pick[v_] := Module[{r},
    r = RandomInteger[range];
    While[MemberQ[excludeZero, v] && r == 0, r = RandomInteger[range]];
    r
  ];
  values = pick /@ vars;
  Thread[vars -> values]
];

EvaluationVector[expr_, rules_] := Module[{val},
  val = Quiet @ Check[Together[expr /. rules], $Failed];
  If[
    val === $Failed ||
      ! FreeQ[val, ComplexInfinity | Indeterminate | DirectedInfinity],
    $Failed,
    Flatten[{val}]
  ]
];

Options[ProbeMatrix] = {
  "Samples" -> Automatic,
  "Variables" -> {J, m1, m2, z, a, b},
  "ExtraAtoms" -> {},
  "AtomHeads" -> {delCoeff},
  "IntegerRange" -> {-7, 7},
  "MaxAttempts" -> Automatic,
  "Seed" -> 12345,
  "Verbose" -> True
};

ProbeMatrix[exprs_List, OptionsPattern[]] := Module[
  {
    n, vars0, atoms, vars, samples, maxAttempts, attempts = 0,
    rows = {}, rules, vals, good, seed, verbose, excludeZero
  },
  n = Length[exprs];
  vars0 = OptionValue["Variables"];
  atoms = DependenceAtoms[
    exprs,
    "ExtraAtoms" -> OptionValue["ExtraAtoms"],
    "AtomHeads" -> OptionValue["AtomHeads"]
  ];
  vars = DeleteDuplicates @ Join[vars0, atoms];
  samples = Replace[OptionValue["Samples"], Automatic :> Max[32, n + 12]];
  maxAttempts = Replace[OptionValue["MaxAttempts"], Automatic :> 5*samples];
  seed = OptionValue["Seed"];
  verbose = OptionValue["Verbose"];
  excludeZero = Select[
    vars,
    MemberQ[{J, m1, m2, z, x, mA, u, v}, #] &
  ];
  If[seed =!= Automatic, SeedRandom[seed]];

  While[Length[rows] < samples && attempts < maxAttempts,
    attempts++;
    rules = SafeRationalRules[
      vars,
      "IntegerRange" -> OptionValue["IntegerRange"],
      "ExcludeZeroFor" -> excludeZero,
      "Seed" -> Automatic
    ];
    vals = EvaluationVector[#, rules] & /@ exprs;
    good = FreeQ[vals, $Failed] &&
      SameQ @@ (Length /@ vals) &&
      AllTrue[Flatten[vals], NumberQ];
    If[
      good,
      rows = Take[Join[rows, Transpose[vals]], UpTo[samples]]
    ];
  ];

  If[verbose,
    Print["Number of constraints: ", n];
    Print["Number of detected delCoeff-like atoms: ", Length[atoms]];
    Print["Number of scalar probe rows: ", Length[rows]];
    Print["Probe attempts used: ", attempts];
  ];

  <|
    "Matrix" -> rows,
    "VariablesUsed" -> vars,
    "DetectedAtoms" -> atoms,
    "Rows" -> Length[rows],
    "Attempts" -> attempts
  |>
];

(* Pivot columns of one row reduction are an independent constraint set. *)
RankIncrementIndependentColumns[m_] := Module[{rr, nonzeroRows},
  rr = RowReduce[m];
  nonzeroRows = Select[rr, AnyTrue[#, Not @* PossibleZeroQ] &];
  First @ FirstPosition[#, _?(Not @* PossibleZeroQ)] & /@ nonzeroRows
];

Options[VerifyRelation] = {
  "SimplifyFunction" -> FullSimplify,
  "Assumptions" -> True
};

VerifyRelation[exprs_List, coeffs_List, OptionsPattern[]] := Module[
  {simp, assump, rel, fast},
  simp = OptionValue["SimplifyFunction"];
  assump = OptionValue["Assumptions"];
  rel = Total[MapThread[#1*#2 &, {coeffs, exprs}]];
  fast = Quiet @ Check[Together[rel], $Failed];
  If[
    fast =!= $Failed && TrueQ[fast === 0],
    True,
    Quiet @ Check[
      TrueQ[simp[If[fast === $Failed, rel, fast] == 0, assump]],
      False
    ]
  ]
];

RelationToRule[coeffs_List, names_List] := Module[{nz},
  nz = Flatten @ Position[coeffs, _?(# =!= 0 &)];
  <|
    "NonzeroIndices" -> nz,
    "NonzeroNames" -> names[[nz]],
    "Coefficients" -> coeffs[[nz]],
    "Formula" -> HoldForm[
      Total[MapThread[#1*#2 &, {coeffs[[nz]], names[[nz]]}]] == 0
    ]
  |>
];

Options[CheckNullDependenceByProbing] = Join[
  Options[ProbeMatrix],
  Options[VerifyRelation],
  {
    "Names" -> Automatic,
    "ProbeExpressions" -> Automatic,
    "VerifySymbolically" -> True
  }
];

CheckNullDependenceByProbing[exprs_List, OptionsPattern[]] := Module[
  {
    names, probeExprs, probe, m, rank, null, independent, dependent,
    verified, relations
  },
  names = Replace[
    OptionValue["Names"],
    Automatic :> Array[Subscript[f, #] &, Length[exprs]]
  ];
  probeExprs = Replace[
    OptionValue["ProbeExpressions"],
    Automatic :> exprs
  ];

  If[
    Length[names] =!= Length[exprs],
    Return[<|"Error" -> "Names and constraints have different lengths."|>]
  ];
  If[
    Length[probeExprs] =!= Length[exprs],
    Return[<|"Error" -> "ProbeExpressions and constraints have different lengths."|>]
  ];

  probe = ProbeMatrix[
    probeExprs,
    "Samples" -> OptionValue["Samples"],
    "Variables" -> OptionValue["Variables"],
    "ExtraAtoms" -> OptionValue["ExtraAtoms"],
    "AtomHeads" -> OptionValue["AtomHeads"],
    "IntegerRange" -> OptionValue["IntegerRange"],
    "MaxAttempts" -> OptionValue["MaxAttempts"],
    "Seed" -> OptionValue["Seed"],
    "Verbose" -> OptionValue["Verbose"]
  ];
  m = probe["Matrix"];
  If[
    Length[m] == 0,
    Return[<|"Error" -> "No valid numeric probe rows were generated."|>]
  ];

  independent = RankIncrementIndependentColumns[m];
  rank = Length[independent];
  null = If[rank == Length[exprs], {}, NullSpace[m]];
  dependent = Complement[Range[Length[exprs]], independent];

  relations = RelationToRule[#, names] & /@ null;
  If[
    TrueQ[OptionValue["VerifySymbolically"]],
    verified = VerifyRelation[
        exprs,
        #,
        "SimplifyFunction" -> OptionValue["SimplifyFunction"],
        "Assumptions" -> OptionValue["Assumptions"]
      ] & /@ null,
    verified = ConstantArray[Missing["NotVerified"], Length[null]]
  ];
  relations = MapThread[
    Append[#1, "VerifiedSymbolically" -> #2] &,
    {relations, verified}
  ];

  <|
    "Summary" -> <|
      "NumberOfConstraints" -> Length[exprs],
      "ProbeRows" -> Length[m],
      "ProbeRank" -> rank,
      "Nullity" -> Length[exprs] - rank,
      "IndependentCount" -> Length[independent],
      "DependentCount" -> Length[dependent],
      "AllCandidateRelationsVerified" ->
        And @@ Replace[verified, {} -> {True}]
    |>,
    "IndependentIndices" -> independent,
    "IndependentNames" -> names[[independent]],
    "DependentIndices" -> dependent,
    "DependentNames" -> names[[dependent]],
    "Relations" -> relations,
    "Probe" -> probe,
    "Matrix" -> m
  |>
];

nullDepReport = CheckNullDependenceByProbing[
  lst,
  "ProbeExpressions" -> (
    lst /. {
      x -> (u^2 + v^2)^2,
      mA -> (u^2 - v^2)/2
    }
  ),
  "Variables" -> {J, m1, m2, u, v, a, b},
  "AtomHeads" -> {delCoeff},
  "IntegerRange" -> {-11, 11},
  "Seed" -> 20260629,
  "VerifySymbolically" -> True,
  "Assumptions" ->
    Element[{J, m1, m2, x, a, b, mA}, Reals] && x > 4*mA^2
];

nullDepReport["Summary"];
nullDepReport["DependentNames"];
nullDepReport["Relations"][[
  All,
  {
    "NonzeroIndices", "NonzeroNames", "Coefficients", "Formula",
    "VerifiedSymbolically"
  }
]];

lstIndependent = lst[[nullDepReport["IndependentIndices"]]];
Length /@ {lst, lstIndependent}

(* ===== Input cell 23 ===== *)
nullDepReport["DependentNames"]



