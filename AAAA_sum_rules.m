(* ::Package:: *)

(* Exported from AAAA_sum_rules.original.nb. *)
(* All original Input cells are included in their original order. *)
(* Only the Linear Dependence block has been replaced by the optimized checker. *)

(* ===== Input cell 1 ===== *)
PartialWaveD[d_,J_,x_]:=Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-x)/2];


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
    -gAAA^2(1/(s-mA^2)+1/(t-mA^2)+1/(u-mA^2))+Total[
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
SumRule[s_,t_,s1_,s2_,k_,J_,mA_]:=((1)/((s-s1)*((s-s1)*(s-s2))^(k/2))-1/((4*mA^2-s-t-s1)*((4*mA^2-s-t-s1)*(4*mA^2-s-t-s2))^(k/2)))*PartialWaveD[10,J,1+(2*t)/(s-4*mA^2)];

(* ===== Input cell 7 ===== *)
ToString[SeriesCoefficient[SumRule[x,t,s1,s2,0,J,mA],{t,0,0}], InputForm, PageWidth->Infinity]



(* ===== Input cell 8 ===== *)
DblCtr[mA_, k_Integer, q_Integer, Nmax_Integer] := -Residue[ Residue[ 1/(s*t) ( Mfwdlow[s, t, mA, Nmax]/(s^(k-q)*t^q) - Mfwdlow[t, s, mA, Nmax]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

(* ===== Input cell 9 ===== *)
(* Double Contour Integral *)
DblCtrSum[s_,J_,k_,q_]:=Residue[((1)/(s*t)*((1)/(s^(k-q)*t^q)-(1)/(s^q*t^(k-q)))-(1)/((4*mA^2-s-t)*t)*((1)/((4*mA^2-s-t)^(k-q)*t^q)-(1)/((4*mA^2-s-t)^q*t^(k-q))))*PartialWaveD[10,J,1+(2*t)/(s-4*mA^2)],{t,0}];

xkq[z_,J_,k_,q_] := DblCtrSum[z,J,k,q];


(* ---------- index list: k, k-q > q ---------- *)
kMax = 13;
nullIndex = Flatten[
  Table[{k, q}, {k, 1, kMax}, {q, 0, Floor[(k - 1)/2]}], 1];
(* literal pseudo code would be {q, 0, Floor[k/2]}: adds 10 identically-zero rows *)
Length[nullIndex] 


MlowPole[s_, t_, mA_] := -1^2 (1/(s - mA^2) + 1/(t - mA^2) + 1/(4mA^2 - s - t - mA^2));
DblCtrFast[k_Integer, q_Integer] := Residue[Residue[
   1/(s t) (MlowPole[s, t, m]/(s^(k - q) t^q) - MlowPole[t, s, m]/(t^(k - q) s^q)),
   {s, Infinity}], {t, 0}];
   
DistributeDefinitions[xkq, DblCtrSum, PartialWaveD, DblCtrFast, MlowPole];


LaunchKernels[];
nullX = ParallelMap[Together@FunctionExpand@xkq[z, J, #[[1]], #[[2]]] &, nullIndex];
nullA = ParallelMap[Together@DblCtrFast[#[[1]], #[[2]]] &, nullIndex];




(* ===== Input cell 19 ===== *)
ClearAll[largeJ];



largeJ[list_List,J_Symbol:J]:=Module[{p},p=Max[Exponent[Together@Simplify[#],J]&/@list];
	Assuming[J>0,Limit[list/J^p,J->Infinity]]];

Jlimit = largeJ[nullX,J]

Length[Jlimit]


ToString[nullA,InputForm,PageWidth->Infinity]


ToString[nullX,InputForm,PageWidth->Infinity]


(* Exact polynomial coefficient rank; earlier labels are kept first. *)
(* Also returns the RREF rows: column d satisfies                      *)
(*   exprs[[d]] == rr[[All, d]] . exprs[[keep]]   (constant coefficients) *)
FindRedundantNullsLabeled[exprs_List, labels_List, vars_List] := Module[
  {ex, den, polys, monos, mat, rr, keep, drop},
  ex    = Together /@ exprs;
  den   = Fold[PolynomialLCM, 1, Denominator /@ ex];
  polys = Expand[Cancel[den #]] & /@ ex;
  monos = Union @@ (First /@ CoefficientRules[#, vars] & /@ polys);
  mat   = Table[Fold[Coefficient[#1, #2[[1]], #2[[2]]] &, p, Transpose[{vars, pw}]],
            {pw, monos}, {p, polys}];
  rr    = Select[RowReduce[mat], AnyTrue[#, # != 0 &] &];
  keep  = First@FirstPosition[#, v_ /; v != 0] & /@ rr;
  drop  = Complement[Range@Length@labels, keep];
  <|"Keep" -> keep, "Drop" -> drop, "Coeff" -> rr,
    "KeepLabels" -> labels[[keep]], "DropLabels" -> labels[[drop]]|>
];

rep = FindRedundantNullsLabeled[nullX, nullIndex, {z, J}];
keepPos = rep["Keep"];  dropPos = rep["Drop"];  cf = rep["Coeff"];
rep["DropLabels"]
