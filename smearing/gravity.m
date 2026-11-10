FileDir = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{FileDir, "..", "SDPB.m"}]];

(* partial-wave basis *)
d = 6;
calP[J_, z_] := Hypergeometric2F1[-J, J + d - 3, (d - 2)/2, (1 - z)/2];

(* improved sum rules *)
kappa[J_] :=
  Together[(2 + v)/(1 + v)^2 calP[J, 1 + 2 v] -
    v^2 ((4 + 3 v)/(1 + v)^2 calP[J, 1] + (4 v)/(
        1 - v^2) (D[calP[J, z], z] /. {z -> 1}))];
(* null constraints *)
chi[k_, J_] :=
  Together[
   Module[{w}, (2 +
        v) calP[J, 1 + 2 v]/(v (1 + v))/(v (1 + v))^(k/2) -
     SeriesCoefficient[(2 + w) (1 -
         w) (1 + 2 w)/((v - w) w (1 - v) (1 + w) (1 + v + w))*
       calP[J, 1 + 2 w]/(w (1 + w))^(k/2), {w, 0, -1}]]];

QR[expr_] := Module[{num = Numerator[expr], den = Denominator[expr], R},
    If[FreeQ[den, v], {CoefficientList[expr, v], 0},
        R = Together[(1 - v) expr] /. v -> 1;
        {CoefficientList[PolynomialQuotient[num - R den/(1 - v), den, v], v], R}
    ]
];

(* ======================================================================
   CACHED SMEARING
   For a kernel K(v) = sum_j q_j v^j + R/(1-v), v = -p^2/m2, y = 1/m2:
     int_0^1 p^n K dp = sum_j M[n,j] y^j + R I_n(y),
     M[n,j] = q_j (-1)^j/(n+2j+1)                     (exact, depends on J only)
     I_n(y) = int_0^1 p^n/(1+y p^2) dp
            = 2F1(1,(n+1)/2;(n+3)/2;-y)/(n+1)          (depends on m2 only)
   Level 1: kernelData[J] -> the three QR divisions are done ONCE per J
            (the old rowF/rowX redid them 7+6+4 = 17 times per state).
   Level 2: Ivals -> 4 hypergeometrics per mass instead of 17, via the
            downward recurrence I_n = 1/(n+1) - y I_{n+2}. This is
            stable for 0 < y <= 1 (m2 >= 1): an error in I_{n+2}
            is multiplied by y <= 1.
   Level 3: the Chebyshev masses are the same for all J, so their I_n
            are computed once and shared by every spin.
   Checked against the old rowF/rowX at 200 digits (J = 0,2,4,6,10,44;
   m2 = 1 ... 1e6): identical to working precision.
   ====================================================================== *)
fList = Range[3/2, 15/2]; iList4 = Range[0, 5]; iList6 = Range[0, 3];

polyMat[{qc_, R_}, nList_] :=
  {Table[qc[[j + 1]] (-1)^j/(n + 2 j + 1), {n, nList}, {j, 0, Length[qc] - 1}], R};

kernelData[J_] := <|
  "F"  -> polyMat[QR[kappa[J]], fList],
  "X4" -> polyMat[QR[chi[4, J]], iList4],
  "X6" -> polyMat[QR[chi[6, J]], iList6]|>;

(* nList must be a run of step 1 (both parities present), as fList and iList4 are *)
Ivals[nList_, y_] := Module[{desc = Reverse[Sort[nList]], tab = <||>},
  Scan[(tab[#] = Hypergeometric2F1[1, (# + 1)/2, (# + 3)/2, -y]/(# + 1)) &, Take[desc, 2]];
  Scan[If[! KeyExistsQ[tab, #], tab[#] = 1/(# + 1) - y tab[# + 2]] &, Drop[desc, 2]];
  tab /@ nList];

(* chi_{k,J} = 0 for J < k gives an empty coefficient matrix: handled explicitly *)
evalBlock[{M_, R_}, y_, Iv_] := Module[
  {dd = If[M === {} || First[M] === {}, 0, Last[Dimensions[M]]]},
  If[dd == 0, ConstantArray[0, Length[Iv]], M . (y^Range[0, dd - 1])] + R Iv];

(* one state: {F (7), X4 (6), X6 (4)} with the mass prefactors m^-4, m^-14, m^-20 *)
stateRowI[kd_, m2_, {IF_, IX_}] := With[{y = 1/m2},
  Join[evalBlock[kd["F"], y, IF]/m2^2,
       evalBlock[kd["X4"], y, IX] m2^(-7),
       evalBlock[kd["X6"], y, Take[IX, Length[iList6]]] m2^(-10)]];
stateRow[kd_, m2_] := stateRowI[kd, m2, {Ivals[fList, 1/m2], Ivals[iList4, 1/m2]}];

stBlock[row_] := PositiveMatrixWithPrefactor[DampedRational[1, {}, 1/E, x], {{row}}];


LaunchKernels[];


g2Fixed = 0;
bmax = 40;
dphi = Pi/4;

(* ---- refined sampling parameters (tune here) ---- *)
hiJmin = 44; hiJmax = 120;
hiBmin = 17; hiBmax = 22; hiBstep = 1/20;
refJmin = 38; refJmax = 42;
dphiFine = Pi/16;
refBmin = 19; refBmax = 20;
fineBmin = 99/5; fineBmax = 20; fineBstep = 1/1000;

(* ---- hooks for adaptive refinement: states/b-points found negative by a
        positivity scan of the previous solution are appended here ---- *)
extraStates = {{6, 3813676789666033/2205434564343156}, {8, 1365884739111897/1232843844564337}, {8, 1582762659416429/626948676181969}, {10, 1039594443994147/619487383610101}, {12, 1420883884224059/610805280875514}, {14, 3701170151192694/618900355491689}, {16, 2561687544364103/341617995739014}, {18, 3992577713709449/433796791226731}, {20, 22772157853128586/4162991224784495}, {20, 1159038419082351/104457428707645}, {22, 1548540799414522/241295207777833}, {24, 2213825107984081/297546535042458}, {120, 3779777465457455/25260924769092}, {130, 13765184671879869/78588649348060}, {150, 20773327532127404/89432730356097}, {200, 20876548008012649/50858742604489}, {300, 19803485912660804/21561886812077}};
extraB = {2656768542969163/133498269154267};
(* written and merged round by round by checkPositivity.m *)
If[FileExistsQ[FileNameJoin[{FileDir, "extraConstraints.m"}]],
  Get[FileNameJoin[{FileDir, "extraConstraints.m"}]];
  Print["loaded extraConstraints.m: ", Length[extraStates], " states, ", Length[extraB], " b points"]];

PMP2SDP[datfile_, prc_: 1000] := Module[
    {
        npts = 150, bBlock, bigB, mmax, nNull,
        largeBBlock, Acoef, trunc, Nexp, ent, osc, nullZeros, cG, c2, c3,
        norm, obj, spins, hiSpins, allSpins, kdList, bList
    },

    chebMasses[np_, prec_] := Module[{phi = Table[Pi/2 + Pi/2 Cos[(k + 1/2) Pi/np], {k, 0, np - 1}]},
        N[1/(1 - Sin[#/2]^2) & /@ phi, prec]
    ];
    suppMasses[J_, dphi_, bmax_, prec_] := Module[{nn = J + 3/2, K},
        K = Floor[ N[Min[nn Pi, bmax]/dphi, 50] ];
        Table[N[1/Sin[j dphi/(2 nn)]^2, prec], {j, 1, K}]
    ];
    refMasses[J_, prec_] := If[refJmin <= J <= refJmax,
        Module[{nn = J + 3/2, kmin, kmax},
            kmin = Ceiling[ N[refBmin/dphiFine, 50] ];
            kmax = Floor[ N[Min[refBmax, nn Pi, bmax]/dphiFine, 50] ];
            Table[N[1/Sin[k dphiFine/(2 nn)]^2, prec],
                {k, Select[Range[kmin, kmax], !IntegerQ[# dphiFine/dphi] &]}]
        ],
        {}
    ];
    hiMasses[J_, prec_] := Table[N[(2 J/b)^2, prec], {b, Range[hiBmin, hiBmax, hiBstep]}];

    bigB = 40; mmax = 2; nNull = 10;

    (* ---- Level 1 cache: exact per-J data, built in parallel, heaviest J first ---- *)
    spins = Range[0, 42, 2]; hiSpins = Range[hiJmin, hiJmax, 2];
    allSpins = Union[spins, hiSpins, extraStates[[All, 1]]];
    kdList = ParallelMap[kernelData, Reverse[allSpins], Method -> "FinestGrained"];
    kdCache = AssociationThread[Reverse[allSpins] -> kdList];

    (* ---- Level 3 cache: Chebyshev masses and their I_n, shared by all spins ---- *)
    chebM = chebMasses[npts, prc];
    chebI = ParallelMap[{Ivals[fList, 1/#], Ivals[iList4, 1/#]} &, chebM];

    DistributeDefinitions[kdCache, chebM, chebI, fList, iList4, iList6,
        polyMat, Ivals, evalBlock, stateRowI, stateRow, stBlock,
        suppMasses, refMasses, hiMasses, dphi, bmax, dphiFine, refJmin, refJmax,
        refBmin, refBmax, hiBmin, hiBmax, hiBstep];

    Acoef[n_] := 2^n Gamma[(n + 1)/2]/Gamma[(3 - n)/2];
    trunc[e_] := Normal[ Series[e bb^(5/2), {bb, Infinity, mmax}] ]/bb^(5/2);
    Nexp = Max[fList] + 1;

    osc[n_, mmax_] := Module[{ser, c1, c2, al, be},
        ser = Expand[ Normal[ Series[ HypergeometricPFQ[{(n + 1)/2}, {2, (n + 3)/2}, -bb^2/4]/(n + 1),
            {bb, Infinity, mmax + 3} ] ] ] /. E^(y_) :> E^Expand[y];

        c1 = Coefficient[ser, E^(I bb)]; c2 = Coefficient[ser, E^(-I bb)];
        al = c1 E^(5 I Pi/4); be = c2 E^(-5 I Pi/4);
        {FullSimplify[ ExpandAll[al + be] ], FullSimplify[ExpandAll[I (al - be)] ] }
    ];

    ent[n_] := Module[{Bn, Cn},
        {Bn, Cn} = trunc /@ osc[n, mmax];
        Expand[bb^Nexp {Acoef[n] bb^(-n - 1) + Bn, Cn, Acoef[n] bb^(-n - 1) - Bn} /. bb -> bigB + x]
    ];

    (* b-blocks evaluated numerically in parallel. In the old code they stayed as exact
       HypergeometricPFQ[...] and were all evaluated serially by N[pols, prc] at the end. *)
    bBlock[b_] := PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, x],
        {{Join[ Table[ N[HypergeometricPFQ[ {(n + 1)/2}, {2, (n + 3)/2}, -b^2/4]/(n + 1), prc], {n, fList} ],
                ConstantArray[0, nNull] ]}}
    ];
    largeBBlock := Module[{e = ent /@ fList, z = ConstantArray[0, nNull]},
        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, x],
            {{Join[e[[All, 1]], z], Join[e[[All, 2]], z]},
            {Join[e[[All, 2]], z], Join[e[[All, 3]], z]}}
        ]
    ];
    bList = Union[Range[1/250, 40, 1/32], Range[fineBmin, fineBmax, fineBstep], extraB];
    DistributeDefinitions[bBlock];   (* prc is already substituted literally; nNull goes along as a dependency *)

    (* Same set of constraints as before (block order differs, which does not matter),
       plus the extraStates / extraB hooks. *)
    pols = Join[
        (* Chebyshev masses: precomputed I_n *)
        Flatten[ParallelTable[stBlock[stateRowI[kdCache[J], chebM[[k]], chebI[[k]]]],
            {J, spins}, {k, Length[chebM]}]],
        (* supplementary + refined masses *)
        Flatten[ParallelTable[stBlock[stateRow[kdCache[J], m2]],
            {J, spins}, {m2, Join[suppMasses[J, dphi, bmax, prc], refMasses[J, prc]]}]],
        (* high spins *)
        Flatten[ParallelTable[stBlock[stateRow[kdCache[J], m2]],
            {J, hiSpins}, {m2, hiMasses[J, prc]}]],
        (* adaptive extra states *)
        ParallelMap[stBlock[stateRow[kdCache[#[[1]]], N[#[[2]], prc]]] &, extraStates],
        ParallelMap[bBlock, bList],
        {largeBBlock}
    ];

    nullZeros = ConstantArray[0, nNull];

    cG = Join[Table[1/(n - 1), {n, fList}], nullZeros];
    c2 = Join[Table[2/(n + 1), {n, fList}], nullZeros];
    c3 = Join[Table[1/(n + 3), {n, fList}], nullZeros];
    norm = -c3;
    obj = -(cG + g2Fixed c2);
    WritePmpJson[datfile, SDP[obj, norm, N[pols, prc] ], prc, getAnalyticSampleData];
];


PMP2SDP["n_pmp.json", 1000]