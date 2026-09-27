(* ::Package:: *)

(* Exported from BBBB_sum_rules.original.nb. *)
(* All original Input cells are included in their original order. *)
(* Only the Linear Dependence block has been replaced by the optimized checker. *)

(* ===== Input cell 1 ===== *)
del[a_, b_] := delCoeff @@ Sort[{a,b}];



validTriples[Nmax_Integer]:=
	Flatten[Table[{a,b}, {a, 0, Nmax}, {b, 0, Nmax-a}], 1];



Mlow[s_, t_, m1_,m2_] := d0+d2*(s^2+t^2+u^2)+d3*(s*t*u)+d4*(s^2+t^2+u^2)^2/.{u -> -s-t};

(* ===== Input cell 2 ===== *)
s1 = 0;


s2 = 0-t;


Ker[s_,t_,m1_,m2_,k_,s1_,s2_] := (Mlow[s,t,m1,m2])/((s-s1)*((s-s1)*(s-s2))^(k/2));

(* ===== Input cell 3 ===== *)
Series[-Residue[Ker[s,t,m1,m2,2,s1,s2],{s,Infinity}],{t,0,2}]

(* ===== Input cell 4 ===== *)
Series[-Residue[Ker[s,t,m1,m2,0,s1,s2],{s,Infinity}],{t,0,4}]//FullSimplify

(* ===== Input cell 5 ===== *)
Series[-Residue[Ker[s,t,m1,m2,4,s1,s2],{s,Infinity}],{t,0,0}]//FullSimplify

(* ===== Input cell 6 ===== *)
PartialWaveD[d_,J_,x_]:=Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-x)/2];


PartialWaveD[4,J,x]-LegendreP[J,x]//FullSimplify

(* ===== Input cell 7 ===== *)
SumRule[s_,t_,s1_,s2_,k_,J_]:=(1/((s-s1)*((s-s1)*(s-s2))^(k/2))-1/((-s-t-s1)*((-s-t-s1)*(-s-t-s2))^(k/2)))*((s)^((10-3)/(2)))/(Sqrt[s])*PartialWaveD[10,J,1+(2*t)/(s)];

(* ===== Input cell 8 ===== *)
Series[SumRule[x,t,s1,s2,0,J],{t,0,0}]

(* ===== Input cell 9 ===== *)
DblCtr[k_Integer, q_Integer] := Residue[ Residue[ 1/(s*t) ( Mlow[s,t,m1,m2]/(s^(k-q)*t^q) -  Mlow[t,s,m1,m2]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

(* ===== Input cell 10 ===== *)
DblCtr[1,0]


DblCtr[2,0]


DblCtr[2,1]


DblCtr[3,0]


DblCtr[3,1]

(* ===== Input cell 11 ===== *)
(* Double Contour Integral *)
DblCtrSum[s_,J_,k_,q_]:=Residue[((1)/(s*t)*((1)/(s^(k-q)*t^q)-(1)/(s^q*t^(k-q)))-(1)/((-s-t)*t)*((1)/((-s-t)^(k-q)*t^q)-(1)/((-s-t)^q*t^(k-q))))*((s)^((10-3)/(2)))/(Sqrt[s])*PartialWaveD[10,J,1+(2*t)/(s)],{t,0}];
xkq[z_,J_,k_,q_] := DblCtrSum[z,J,k,q];

(* ===== Input cell 12 ===== *)
x10 = xkq[z,J,1,0]//FullSimplify


x20 = xkq[z,J,2,0]//FullSimplify


x30 =  xkq[z,J,3,0]//FullSimplify


x40 =  xkq[z,J,4,0]//FullSimplify


x41 = xkq[z,J,4,1]//FullSimplify

(* ===== Input cell 14 ===== *)
x50 = xkq[z,J,5,0]//FullSimplify


x51 = xkq[z,J,5,1]//FullSimplify

(* ===== Input cell 15 ===== *)

x60 = xkq[z,J,6,0]//FullSimplify
x61 = xkq[z,J,6,1]//FullSimplify

(* ===== Input cell 16 ===== *)
x70 = xkq[z,J,7,0]//FullSimplify


x71 = xkq[z,J,7,1]//FullSimplify


x73 = xkq[z,J,7,3]//FullSimplify

(* ===== Input cell 17 ===== *)
x80 = xkq[z,J,8,0]//FullSimplify


x81 = xkq[z,J,8,1]//FullSimplify


x83 = xkq[z,J,8,3]//FullSimplify

(* ===== Input cell 18 ===== *)
x90 = xkq[z,J,9,0]//FullSimplify


x91 = xkq[z,J,9,1]//FullSimplify


x93 = xkq[z,J,9,3]//FullSimplify

(* ===== Input cell 19 ===== *)
x100 = xkq[z,J,10,0]//FullSimplify


x101 = xkq[z,J,10,1]//FullSimplify


x103 = xkq[z,J,10,3]//FullSimplify


x104 = xkq[z,J,10,4]//FullSimplify

(* ===== Input cell 20 ===== *)
x110 = xkq[z,J,11,0]//FullSimplify


x111 = xkq[z,J,11,1]//FullSimplify


x113 = xkq[z,J,11,3]//FullSimplify


x114 = xkq[z,J,11,4]//FullSimplify

(* ===== Input cell 21 ===== *)
x120 = xkq[z,J ,12,0]//FullSimplify


x121 = xkq[z,J ,12,1]//FullSimplify


x123 = xkq[z,J ,12,3]//FullSimplify


x124 = xkq[z,J ,12,4]//FullSimplify

(* ===== Input cell 22 ===== *)
x130 = xkq[z, J,13,0]//FullSimplify


x131 = xkq[z, J,13,1]//FullSimplify


x132 = xkq[z,J,13,2]//FullSimplify


x133 = xkq[z, J,13,3]//FullSimplify


x134 = xkq[z, J,13,4]//FullSimplify

(* ===== Input cell 23 ===== *)
x140 = xkq[z, J,14,0]//FullSimplify


x141 = xkq[z, J,14,1]//FullSimplify


x142 = xkq[z, J,14,2]//FullSimplify


x143 = xkq[z, J,14,3]//FullSimplify


x144 = xkq[z, J,14,4]//FullSimplify

(* ===== Input cell 24 ===== *)
x150 = xkq[z,J,15,0]//FullSimplify;


x151 = xkq[z, J,15,1]//FullSimplify;


x152 = xkq[z, J,15,2]//FullSimplify;


x153 = xkq[z, J,15,3]//FullSimplify;


x154 = xkq[z, J,15,4]//FullSimplify;

(* ===== Input cell 25 ===== *)
x160 = xkq[z,J,16,0]//FullSimplify;


x161 = xkq[z,J,16,1]//FullSimplify;


x162 = xkq[z,J,16,2]//FullSimplify;


x163 = xkq[z,J,16,3]//FullSimplify;


x164 = xkq[z,J,16,4]//FullSimplify;


x165 = xkq[z,J,16,5]//FullSimplify;

(* ===== Input cell 26 ===== *)
x170 = xkq[z,J,17,0]//FullSimplify;

(* ===== Input cell 27 ===== *)
ClearAll[largeJ];



largeJ[list_List,J_Symbol:J]:=Module[{p},p=Max[Exponent[Together@Simplify[#],J]&/@list];
Assuming[J>0,Limit[list/J^p,J->Infinity]]]

(* ===== Input cell 28 ===== *)
pref = z^7;



lst = {x10,x20,x30, x40, x41 ,x50, x51,x60,x61,x70,x71,x73,x80,x81,x83,x90,x91,x93,x100,x101,x103,x104,x110,x111,x113,x114,x120,x121,x123,x124,x130};



Length[lst]



largeJ[lst]//FullSimplify

(* ===== Input cell 29 ===== *)
lst//FullSimplify

(* ===== Input cell 30 ===== *)
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
  "Variables" -> {J, m1, m2, z, a, b},
  "AtomHeads" -> {delCoeff},
  "IntegerRange" -> {-11, 11},
  "Seed" -> 20260629,
  "VerifySymbolically" -> True,
  "Assumptions" -> Element[{J, m1, m2, z, a, b}, Reals]
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

(* ===== Input cell 31 ===== *)
nullDepReport["DependentNames"]





