(* ::Package:: *)

(* Exported from BBBB_sum_rules.original.nb. *)
(* All original Input cells are included in their original order. *)
(* Only the Linear Dependence block has been replaced by the optimized checker. *)



del[a_, b_] := delCoeff @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];

(* ---------- amplitude ansatz ---------- *)

(* Commented-out single-term prototype kept for reference *)

(* BBBB scattering change the ansatz to be s-u symmetric *)
Mlow[s_, t_, mA_, Nmax_Integer] := Module[{u},
	u = -s-t;
    gAB^2(1/(s-mA^2)+1/(t-mA^2)+1/(u-mA^2))+Total[
        Function[{ab},
            del[ab[[1]], ab[[2]]]
            * (s)^ab[[1]]
            * (u)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] ];



(* ===== Input cell 2 ===== *)
s1 = 0;


s2 = 0-t;


Ker[s_,t_,mA_,k_,s1_,s2_,Nmax_Integer] := (Mlow[s,t,mA,Nmax])/((s-s1)*((s-s1)*(s-s2))^(k/2));

(* 0 subtraction *)
SeriesCoefficient[-Residue[Ker[s,t,mA,0,s1,s2,10],{s,Infinity}],{t,0,0}]


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
x00 = xkq[z,J,0,0]//FullSimplify
x10 = xkq[z,J,1,0]//FullSimplify
x11 = xkq[z,J,1,1]//FullSimplify


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



lst = {x00,x10,x11,x20,x30, x40, x41 ,x50, x51,x60,x61,x70,x71,x73,x80,x81,x83,x90,x91,x93,x100,x101,x103,x104,x110,x111,x113,x114,x120,x121,x123,x124,x130};



Length[lst]



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
  {x00,x10,x11,x20,x30,x40,x41,x50,x51,x60,x61,x70,x71,x73,
   x80,x81,x83,x90,x91,x93,x100,x101,x103,x104,x110,x111,
   x113,x114,x120,x121,x123,x124,x130},
  {z, J}
];

nullDepReport

