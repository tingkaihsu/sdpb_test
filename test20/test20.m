(* ::Package:: *)

test20Directory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{test20Directory, "..", "SDPB.m"}]];

(* ------------------------------------------------------------------------ *)
(* test20: mixed AA / AB / BB dual bootstrap, objective X_H.                 *)
(*                                                                          *)
(*   A : Z2-even scalar, mass mA (external, and first isolated state J = 0)  *)
(*   B : Z2-odd massless scalar (external, and a light state in AB, J = 0)   *)
(*   H : neutral spin-2 isolated state at m2sq = 1                           *)
(*   continuum : x >= mgap; neutral (even J) and AB (all J)                  *)
(*                                                                          *)
(*   test16 idea: the light states A and B are ordinary states of the       *)
(*   spectrum; same kernels, positive weights gAAA^2, gABB^2, gAAA gABB.     *)
(*                                                                          *)
(*   Objective : X_H = gHAA gHBB (off-diagonal of the H block), in RESIDUE *)
(*               units (option R): the X_H marker carries phiX like the    *)
(*               rest of its element.  Density units: X_H/phiX[m2sq].      *)
(*                                                                          *)
(*   Heavy states (H, continuum) carry the D=10 phase-space weights on     *)
(*   every component; the light block A/B carries weight 1 (residues).    *)
(*   Normalization : g0AA = rAA, g0BB = 1 (all states, A included)          *)
(*   Bound : xhSign X_H / g0BB <= alpha_AA rAA + alpha_BB                    *)
(*   Invariant reported : X_H / Sqrt[g0AA g0BB] = bound / Sqrt[rAA]          *)
(*                                                                          *)
(*   includeABB = False : gABB = 0 (symmetric control; the bound / Sqrt[rAA] *)
(*                        must not depend on rAA)                           *)
(*   includeABB = True  : gABB free (A in BB -> BB, B in AB -> AB, A in     *)
(*                        AA -> BB all switched on)                         *)
(*                                                                          *)
(* Functional coordinates, in order:                                       *)
(*   {g0AA, g0BB, XH, cross nulls (36), AA nulls (31), BB nulls (31)}       *)
(*                                                                          *)
(* Basis of the neutral 2 x 2 blocks: {AA, BB}.  The AB channel is a        *)
(* separate 1 x 1 block for every J.                                        *)
(* ------------------------------------------------------------------------ *)

mAsq = 1/5;
mA = Sqrt[mAsq];
JA = 0;
m2sq = 1;
J2 = 2;
mgap = 2;

rAA = 1;            (* g0AA/g0BB; scan it *)
xhSign = 1;         (* +1: maximize X_H,  -1: maximize -X_H *)
includeABB = False; (* see header *)

crossNullIndices = Range[0, 35];
aaNullIndices = Range[0, 30];
bbNullIndices = Range[0, 30];

nCross = Length[crossNullIndices];
nAA = Length[aaNullIndices];
nBB = Length[bbNullIndices];

list0 = ConstantArray[0, nCross + nAA + nBB];
zC = ConstantArray[0, nCross];
zA = ConstantArray[0, nAA];
zB = ConstantArray[0, nBB];


(* ------------------------------------------------------------------------ *)
(* NULL CONSTRAINTS: TO BE FILLED                                           *)
(*                                                                          *)
(* COUPLING NORMALIZATION.  K[n, x, J] is the contribution of ONE state of  *)
(* mass^2 x and spin J whose amplitude near the pole is                     *)
(*     - g^2 P_J(cos theta)/(s - x)        (physical sign, P_J(1) = 1)      *)
(* i.e. the sum-rule files' DblCtrSum[x, J, k, q] divided by the phase-     *)
(* space factor.  The same K is used for EVERY state.                       *)
(* Heavy blocks multiply each element ij by phi_ij (see PHASE-SPACE WEIGHTS *)
(* below), which reproduces the sum-rule files exactly (DblCtrSum,          *)
(* KerAABB1 with phiAABB, NABAB with wAB).  The light block keeps weight 1: *)
(* the light states sit below threshold, where phi is complex/singular.     *)
(*                                                                          *)
(* Expansion points: AAAA and BBBB at (s, t) = (0, 0); cross at s0Cross.    *)
(* Low-energy amplitudes are pure polynomials: NO IR poles in the ansatz.   *)
(* ------------------------------------------------------------------------ *)

kernelsFilled = False;

s0Cross = (2 mAsq - mgap)/2;   (* cross expansion point; must be < 0 *)

KAAAA[n_, x_, J_] := 0;   (* AA -> AA null constraint n, AA entry          *)
KBBBB[n_, x_, J_] := 0;   (* BB -> BB null constraint n, BB entry          *)
KAABB[n_, x_, J_] := 0;   (* cross null constraint n, neutral off-diagonal *)
KABAB[n_, x_, J_] := 0;   (* cross null constraint n, AB-channel entry     *)


(* ------------------------------------------------------------------------ *)
(* Positive moments: k = 0 forward sum rules (t = 0), test16 convention.    *)
(*                                                                          *)
(*   g0AA = M_AA(s = 0, t = 0)   (u = 4 mA^2; the AAAA expansion point)     *)
(*   g0BB = M_BB(s = 0, t = 0)                                               *)
(*                                                                          *)
(* = sums over ALL states, light A included, with the weights below.       *)
(* Both weights are positive on every state:                                *)
(*   g0AAWeight: heavy x >= m2sq > 4 mA^2;  A: 1/mA^2 - 1/(3 mA^2) > 0      *)
(*   g0BBWeight: 2/x > 0;                   A: 2/mA^2                       *)
(* (The symmetric point s = 2 mA^2 would give the A state -2/mA^2.)         *)
(* B is not exchanged in AA -> AA or BB -> BB (Z2), so it has no g0 entry.  *)
(* ------------------------------------------------------------------------ *)

(* D=10, AA -> AA at s = 0, t = 0, coupling normalization *)
g0AAWeight[x_] := 1/x + 1/(x - 4 mA^2);

(* D=10, BB -> BB at s = 0, t = 0, coupling normalization *)
g0BBWeight[x_] := 2/x;


(* ------------------------------------------------------------------------ *)
(* PHASE-SPACE WEIGHTS (heavy states only), D = 10: phi = (2p)^7/Sqrt[x]    *)
(*                                                                          *)
(* Element ij of a heavy block carries phi_ij on ALL its components         *)
(* (g0AA, g0BB, XH and every null row).  The off-diagonal weight is the     *)
(* geometric mean, so the PSD cone is unchanged (positive diagonal          *)
(* congruence D = diag(Sqrt[phiAA], Sqrt[phiBB])).                          *)
(*   phiBB g0BBWeight = 2 x^2  (the original BB forward weight)             *)
(*   phiX            = phiAABB of ABAB_check.m                              *)
(*   phiAB           = wAB of ABAB_check.m                                  *)
(* ------------------------------------------------------------------------ *)

phiAA[x_] := (x - 4 mA^2)^(7/2)/Sqrt[x];
phiBB[x_] := x^3;
phiX[x_]  := Sqrt[phiAA[x] phiBB[x]];     (* = x^(5/4) (x - 4 mA^2)^(7/4) *)
phiAB[x_] := (x - mA^2)^7/x^4;

(* Heavy neutral state at (x, J): 2 x 2 in {AA, BB}.
   xh = xhSign for H (objective marker), 0 for the continuum.
   Option R: the marker sits INSIDE the phiX-weighted element, so
   Tr[rhoTilde . M] gives XH = xh phiX rhoTildeX = xh gHAA gHBB (residue). *)
neutralEntries[j_, x_, xh_] := Module[{aa, ab, bb},
  aa = phiAA[x] entry[{g0AAWeight[x], 0, 0}, zC, aaRow[KAAAA, x, j], zB];
  ab = phiX[x]  entry[{0, 0, xh/2}, -1/2 crossRow[KAABB, x, j], zA, zB];
  bb = phiBB[x] entry[{0, g0BBWeight[x], 0}, zC, zA, bbRow[KBBBB, x, j]];
  twoByTwo[aa, ab, bb]
];


(* ------------------------------------------------------------------------ *)
(* Helpers                                                                  *)
(* ------------------------------------------------------------------------ *)

(* one matrix entry = list over all functional coordinates *)
entry[head_List, cross_List, aa_List, bb_List] := Join[head, cross, aa, bb];

(* symmetric 2 x 2 matrix of coordinate vectors in the basis {AA, BB} *)
twoByTwo[aa_, ab_, bb_] := {{aa, ab}, {ab, bb}};

crossRow[f_, x_, J_] := Table[f[n, x, J], {n, crossNullIndices}];
aaRow[f_, x_, J_] := Table[f[n, x, J], {n, aaNullIndices}];
bbRow[f_, x_, J_] := Table[f[n, x, J], {n, bbNullIndices}];


LaunchKernels[];


PMP2SDP[datfile_, prec_:600] := Module[
    {
        npts, phiSamples, massSamples, Jmax,
        evenSpinSamples, allSpinSamples,
        PolyA, PolyH, PolyN, PolyAB, pols, norm, obj,
        functionalCount, expectedBlocks
    },
    If[! TrueQ[4 mAsq < Min[m2sq, mgap] && Abs[4 mAsq - m2sq] > 1/20],
      Print["Invalid masses: need 4 mA^2 < Min[m2sq, mgap] and 4 mA^2 away from m2sq."];
      Abort[]
    ];
    (* B (x = 0) must lie strictly above the cross expansion point. *)
    If[! TrueQ[2 mAsq - mgap < s0Cross <= -mAsq^2/(mgap - 2 mAsq)],
      Print["Invalid s0Cross: need 2 mA^2 - mgap < s0Cross <= -mA^4/(mgap - 2 mA^2) < 0."];
      Abort[]
    ];
    (* Light entries must be finite (no light mass at an expansion point). *)
    If[! FreeQ[
        {aaRow[KAAAA, mAsq, JA], bbRow[KBBBB, mAsq, JA],
         crossRow[KAABB, mAsq, JA], crossRow[KABAB, 0, 0],
         g0AAWeight[mAsq], g0BBWeight[mAsq]},
        ComplexInfinity | Indeterminate | DirectedInfinity],
      Print["A light-state kernel is singular. Check the expansion points."];
      Abort[]
    ];
    (* Phase-space weights must be real and positive on every heavy state
       (H at m2sq, continuum x >= mgap; all are increasing for x > 4 mA^2). *)
    If[! AllTrue[
        {phiAA[m2sq], phiBB[m2sq], phiX[m2sq],
         phiAA[mgap], phiBB[mgap], phiX[mgap], phiAB[mgap]},
        TrueQ[Positive[N[#]]] &],
      Print["A phase-space weight is not real and positive on a heavy state."];
      Abort[]
    ];
    If[! TrueQ[kernelsFilled],
      Print["WARNING: null-constraint kernels are placeholders (= 0). ",
        "The PMP shows the structure only; do not run SDPB on it."]
    ];

    (* Paper eq. (D.1): Chebyshev nodes in the conformal angle.
       Here m = 1-mgap/z = Sin[phi/2]^2. *)
    npts = 150;
    phiSamples = N[Table[
      Pi/2 + Pi/2 Cos[(k + 1/2) Pi/npts],
      {k, 0, npts - 1}
    ], prec];
    massSamples = Sin[#/2]^2 & /@ phiSamples;

    If[
      Length[massSamples] =!= npts ||
        ! AllTrue[massSamples, 0 < # < 1 &],
      Print["Invalid mass-sampling grid."];
      Abort[]
    ];

    (* Neutral AA/BB states have even spin; AB states have all spins. *)
    Jmax = 100;
    evenSpinSamples = Range[0, Jmax, 2];
    allSpinSamples = Range[0, Jmax];

    (* Light states, test16 style: A (x = mAsq, J = 0) and B (x = 0, J = 0)
       enter exactly like any other state, with the SAME kernels evaluated
       at their masses and spins.  The low-energy amplitude is a pure
       crossing-symmetric polynomial (no IR poles anywhere).
       The two cubic couplings tie the light exchanges into ONE 2 x 2 block
       in (gAAA, gABB):
         AA : A in AA -> AA                               (gAAA^2)
         BB : A in BB -> BB  +  B in AB -> AB             (gABB^2)
         AB : A in AA -> BB (= A in the AB -> AB t-channel) (gAAA gABB)
       includeABB = False -> gABB = 0, only the AA entry survives. *)
    PolyA[y_] := Module[{aa, ab, bb},
      aa = entry[{g0AAWeight[mAsq], 0, 0}, zC, aaRow[KAAAA, mAsq, JA], zB];
      ab = entry[{0, 0, 0}, -1/2 crossRow[KAABB, mAsq, JA], zA, zB];
      bb = entry[{0, g0BBWeight[mAsq], 0},
             crossRow[KABAB, 0, 0],                  (* B in AB -> AB *)
             zA,
             bbRow[KBBBB, mAsq, JA]];                (* A in BB -> BB *)

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        If[TrueQ[includeABB], twoByTwo[aa, ab, bb], {{aa}}]
      ]
    ];

    (* Second isolated state H: neutral 2 x 2 block at x = m2sq, J = 2,
       phase-space weighted.  The objective X_H is its off-diagonal
       (option R: residue units, X_H = gHAA gHBB). *)
    PolyH[j_, x_, y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      neutralEntries[j, x, xhSign]
    ];

    (* Neutral continuum: 2 x 2 {AA, BB}, even J, phase-space weighted. *)
    PolyN[j_, x_, y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      neutralEntries[j, x, 0]
    ];

    (* AB continuum: 1 x 1, all J, weight phiAB.  Only the cross nulls. *)
    PolyAB[j_, x_, y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      {{phiAB[x] entry[{0, 0, 0}, crossRow[KABAB, x, j], zA, zB]}}
    ];

    (* Make the global helpers available on the parallel kernels. *)
    DistributeDefinitions[
      phiAA, phiBB, phiX, phiAB, neutralEntries, g0AAWeight, g0BBWeight,
      entry, twoByTwo, crossRow, aaRow, bbRow,
      KAAAA, KBBBB, KAABB, KABAB, zC, zA, zB, mA, xhSign,
      crossNullIndices, aaNullIndices, bbNullIndices
    ];

    pols = Flatten[{
      {N[PolyA[x], prec]},
      Flatten[ N[ ParallelTable[ PolyH[i, m2sq, x], {i, J2, J2, 2}], prec] ],
      Flatten[ N[ ParallelTable[ PolyN[i, mgap/(1-m), x], {i, evenSpinSamples}, {m, massSamples}], prec] ],
      Flatten[ N[ ParallelTable[ PolyAB[i, mgap/(1-m), x], {i, allSpinSamples}, {m, massSamples}], prec] ]
    }, 1];

    expectedBlocks = 2 + npts (
      Length[evenSpinSamples] + Length[allSpinSamples]
    );
    If[Length[pols] =!= expectedBlocks,
      Print[
        "Unexpected block count: built ", Length[pols],
        ", expected ", expectedBlocks, "."
      ];
      Abort[]
    ];

    If[! FreeQ[pols, _Complex],
      Print["Complex entries in the PMP. Aborting."];
      Abort[]
    ];

    Print[
      "Built ", Length[pols], " blocks: light block A/B (",
      If[TrueQ[includeABB], "2x2", "1x1"],
      "), H (2x2), neutral even J = 0..", Jmax, " (2x2), AB J = 0..", Jmax, " (1x1)."
    ];

    (* Coordinates {g0AA, g0BB, XH, nulls}.
       norm : the quantity being bounded (XH).
       obj  : g0AA = rAA, g0BB = 1.
       Bound: xhSign X_H / g0BB <= alpha_AA rAA + alpha_BB. *)
    norm = -1 * N[Flatten[{{0, 0, 1}, list0}], prec];
    obj = -1 * N[Flatten[{{rAA, 1, 0}, list0}], prec];

    functionalCount = 3 + Length[list0];
    If[Length[norm] =!= functionalCount || Length[obj] =!= functionalCount,
      Print["Objective/normalization dimension mismatch."];
      Abort[]
    ];

    Print["size of norm = ", Length[norm]];
    Print["size of obj = ", Length[obj]];
    Print["Writing ", datfile, "..."];
    WritePmpJson[datfile, SDP[obj, norm, pols], prec, getAnalyticSampleData];

    Print["Wrote ", datfile, "."]
];

PMP2SDP[FileNameJoin[{test20Directory, "n_pmp.json"}], 1000];
