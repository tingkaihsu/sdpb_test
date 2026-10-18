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
    -1^2(1/(s-mA^2)+1/(t-mA^2)+1/(u-mA^2))+Total[
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
SumRule[s_,t_,s1_,s2_,k_,J_]:=(1/((s-s1)*((s-s1)*(s-s2))^(k/2))-1/((-s-t-s1)*((-s-t-s1)*(-s-t-s2))^(k/2)))*PartialWaveD[10,J,1+(2*t)/(s)];

(* ===== Input cell 9 ===== *)
DblCtr[k_Integer, q_Integer] := Residue[ Residue[ 1/(s*t) ( Mlow[s,t,m,10]/(s^(k-q)*t^q) -  Mlow[t,s,m,10]/(t^(k-q)*s^q) ), {s, Infinity}], {t, 0}];

(* Double Contour Integral *)
DblCtrSum[s_,J_,k_,q_]:=Residue[((1)/(s*t)*((1)/(s^(k-q)*t^q)-(1)/(s^q*t^(k-q)))-(1)/((-s-t)*t)*((1)/((-s-t)^(k-q)*t^q)-(1)/((-s-t)^q*t^(k-q))))*((s)^((10-3)/(2)))/(Sqrt[s])*PartialWaveD[10,J,1+(2*t)/(s)],{t,0}];

xkq[z_,J_,k_,q_] := DblCtrSum[z,J,k,q];


(* ---------- index list: k = 1..20, k-q > q ---------- *)
kMax = 16;
nullIndex = Flatten[
  Table[{k, q}, {k, 1, kMax}, {q, 0, Floor[(k - 1)/2]}], 1];
(* literal pseudo code would be {q, 0, Floor[k/2]}: adds 10 identically-zero rows *)
Length[nullIndex]   (* 110 *)

(* ---------- faster DblCtr: polynomial part cancels identically ---------- *)
MlowPole[s_, t_, mA_] := -1^2 (1/(s - mA^2) + 1/(t - mA^2) + 1/(-s - t - mA^2));
DblCtrFast[k_Integer, q_Integer] := Residue[Residue[
   1/(s t) (MlowPole[s, t, m]/(s^(k - q) t^q) - MlowPole[t, s, m]/(t^(k - q) s^q)),
   {s, Infinity}], {t, 0}];

(* check against the original DblCtr on small k *)
And @@ Table[PossibleZeroQ[DblCtr @@ kq - DblCtrFast @@ kq],
  {kq, Select[nullIndex, #[[1]] <= 6 &]}]

(* ---------- the double loop ---------- *)
LaunchKernels[];
DistributeDefinitions[xkq, DblCtrSum, PartialWaveD, DblCtrFast, MlowPole];

nullX = ParallelMap[Together@FunctionExpand@xkq[z, J, #[[1]], #[[2]]] &, nullIndex];
nullA = ParallelMap[Together@DblCtrFast[#[[1]], #[[2]]] &, nullIndex];


largeJ[list_List,J_Symbol:J]:=Module[{p},p=Max[Exponent[Together@Simplify[#],J]&/@list];
Assuming[J>0,Limit[list/J^p,J->Infinity]]]



Jlimit = largeJ[nullX,J]



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




(* ---------- sanity checks on the index sets ---------- *)
dropExpected = {{3,1},{5,2},{6,2},{7,3},{8,3},{9,3},{9,4},{10,4},{11,4},{11,5},{12,4},{12,5},{13,5},{13,6},{14,5},{14,6},{15,5},{15,6},{15,7},{16,6},{16,7}};
rep["DropLabels"] === dropExpected                          (* True *)
Sort@Join[keepPos, dropPos] === Range@Length@nullIndex      (* True *)
Length[keepPos]                                             (* 63 *)

(* ---------- the same relations must hold on all three lists ---------- *)
relHoldsQ[list_] := And @@ Table[
  PossibleZeroQ@Together[list[[d]] - cf[[All, d]] . list[[keepPos]]], {d, dropPos}];
{relHoldsQ[nullX], relHoldsQ[nullA], relHoldsQ[Jlimit]}     (* {True, True, True} *)

(* ---------- final lists without redundancies (aligned by nullIndexInd) ---------- *)
nullIndexInd = nullIndex[[keepPos]];
nullXind     = nullX[[keepPos]];
nullAind     = nullA[[keepPos]];
JlimitInd    = Jlimit[[keepPos]];

(* the common power J^p in largeJ must not change after the reduction *)
And @@ (PossibleZeroQ /@ Together[JlimitInd - largeJ[nullXind, J]])   (* True *)


ToString[JlimitInd,InputForm, PageWidth->Infinity]
