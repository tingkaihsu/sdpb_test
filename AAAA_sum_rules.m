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
Mfwdlow[s_, t_, mA_, Nmax_Integer] := Module[{u},
	u = 4mA^2-s-t;
    gAA^2(1/(s-mA^2)+1/(t-mA^2)+1/(u-mA^2))+Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s-2*mA^2)^ab[[1]]
            * (u-2*mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] ];

(* ===== Input cell 4 ===== *)
Ker[s_,t_,s1_,s2_,k_, mA_]:=(Mfwdlow[s,t,mA,10])/((s-s1)*((s-s1)*(s-s2))^(k/2));


s1 = 2*mA^2;


s2 = 2*mA^2-t;

(* ===== Input cell 5 ===== *)
(* g0 + ... *)
Series[-Residue[Ker[s,t,s1,s2,0,mA],{s,Infinity}],{t,0,0}]//FullSimplify


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

x00 = DblCtrSum[x,0,0]//FullSimplify;

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
lst = {x00,x10,x20,x30, x40, x41 ,x50, x51,x52,x60,x61,x62,x70,x71,x73,x80,x81,x83,x90,x91,x93,x94,x100,x101,x103,x104,x110,x111,x112,x113,x114,x115};



largeJ[lst]//FullSimplify



Length[lst]

(* ===== Input cell 21 ===== *)
lst//FullSimplify

(* Exact polynomial coefficient rank; earlier labels are kept first. *)
ClearAll[FindRedundantNulls];
SetAttributes[FindRedundantNulls, HoldAll];

FindRedundantNulls[names_List, variables_List] := Module[
  {held, labels, expressions, denominator, polynomials,
   monomials, matrix, reduced, keep, drop},
  held = HoldComplete[names];
  expressions = Together /@ names;
  labels = Table[
    Extract[held, {1, i}, HoldForm],
    {i, Length[expressions]}
  ];

  (* One common denominator preserves constant-coefficient dependence. *)
  denominator = Fold[
    PolynomialLCM, 1, Denominator /@ expressions
  ];
  polynomials = Expand[Cancel[denominator #]] & /@ expressions;
  monomials = Union[
    Flatten[(First /@ CoefficientRules[#, variables]) & /@ polynomials, 1]
  ];

  keep = If[monomials === {}, {},
    matrix = Table[
      Fold[
        Coefficient[#1, #2[[1]], #2[[2]]] &,
        polynomial,
        Transpose[{variables, powers}]
      ],
      {powers, monomials}, {polynomial, polynomials}
    ];
    reduced = Select[RowReduce[matrix], AnyTrue[#, # != 0 &] &];
    (First @ FirstPosition[#, value_ /; value != 0]) & /@ reduced
  ];
  drop = Complement[Range[Length[labels]], keep];

  <|
    "KeepLabels" -> labels[[keep]],
    "DropLabels" -> labels[[drop]]
  |>
];

nullDepReport = FindRedundantNulls[
  {x00, x10,x20,x30, x40, x41 ,x50, x51,x52,x60,x61,x62,x70,x71,x73,x80,x81,x83,x90,x91,x93,x94,x100,x101,x103,x104,x110,x111,x112,x113,x114,x115},
  {x, J,mA}
];

nullDepReport

