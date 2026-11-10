(* ::Package:: *)

(* =====================================================================
   checkPositivity.m

   Reads the SDPB solution vector z from a log (the out.txt + z.txt dump
   used for 1st.log / 2nd.log), checks it, and tests positivity of the
   functional BETWEEN grid points (tutorial Sec. 7.5):

     (A) finite states (m2, J): value m^4 * [ F . a + X4 . b4 + X6 . b6 ]
         on a scan that is uniform in the phase (J+3/2) theta,
         sin(theta/2) = 1/m, plus a logarithmic grid in m2;
     (B) impact parameter: fhat(b) = sum_n a_n 1F2(...)/(n+1), for
         0 < b <= bScanMax (also beyond the analytic block at b = 40).

   Every local minimum of the scan that is close to zero is refined by a
   golden-section search. Minima below -violTol are written as exact
   rationals to extraConstraints.m, MERGED with what is already there, in
   the format read by the extraStates / extraB hooks of gravity_cached.m.

   Usage:
     wolframscript -file checkPositivity.m 2nd.log
   or, in a notebook, set logFile below and evaluate the whole file.

   The kernel definitions below must stay identical to the generator
   (fList, iList4, iList6, kappa, chi, QR, mass prefactors).
   ===================================================================== *)

FileDir = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];

(* ---------------------------- parameters ---------------------------- *)
logFile    = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "2nd.log"];
g2Fixed    = 0;                                   (* must equal the value used to write the PMP *)

fList  = Range[3/2, 15/2];                        (* f(p) powers          *)
iList4 = Range[0, 5];                             (* p^i X_4, i = 0..5    *)
iList6 = Range[0, 3];                             (* p^i X_6, i = 0..3    *)

scanSpins  = Join[Range[0, 120, 2], {130, 150, 200, 300}];
dphiScan   = Pi/64;   (* phase step of the scan; keep <= 1/4 of the finest phase spacing of the PMP grid *)
logM2Max   = 8;       (* extra grid m2 = 10^(0 .. logM2Max) for the large-mass region *)
logStep    = 1/20;
scanPrec   = 150;     (* working precision of the state scan *)

bStep      = 1/50;    (* fhat scan *)
bScanMax   = 200;
bPrec      = 60;

polishBelow = 10^-2;  (* refine every grid local minimum whose value is below this *)
violTol     = 10^-15; (* refined minimum counts as a violation if value < -violTol *)
goldenIters = 80;

extraFile  = FileNameJoin[{FileDir, "extraConstraints.m"}];
reportFile = FileNameJoin[{FileDir, "positivity_scan.csv"}];


(* ===================== Step 1: read z from the log ===================== *)
logPath = If[FileExistsQ[logFile], logFile, FileNameJoin[{FileDir, logFile}]];
If[! FileExistsQ[logPath], Print["Log file not found: ", logFile]; Abort[]];
lines = StringTrim /@ ReadList[logPath, String];     (* StringTrim also removes \r of CRLF files *)

(* "1.23e-45" -> "1.23*^-45"; long decimal literals become arbitrary-precision numbers
   whose precision equals the number of digits written *)
toNum[s_String] := Module[{t, x},
  t = StringReplace[StringTrim[s], RegularExpression["[eE]([+-]?[0-9]+)$"] -> "*^$1"];
  x = Quiet@ToExpression[t];
  If[! NumericQ[x], Print["Cannot parse number: ", StringTake[s, UpTo[60]]]; Abort[]];
  If[Precision[x] === MachinePrecision, SetPrecision[x, 30], x]];

nExpected = Length[fList] + Length[iList4] + Length[iList6];

(* header line of an SDPB vector file: "<N> 1" *)
hdrIdx = Flatten@Position[StringMatchQ[lines, RegularExpression["[0-9]+\\s+1"]], True];
hdrIdx = Select[hdrIdx, ToExpression[First@StringSplit[lines[[#]]]] == nExpected &];
Which[
  hdrIdx === {}, Print["No vector of dimension ", nExpected, " found in ", logPath,
                       ". (y.txt has dimension N-1; this script needs the full functional z.)"]; Abort[],
  Length[hdrIdx] > 1, Print["Several vectors of dimension ", nExpected, " in ", logPath,
                            " (lines ", hdrIdx, "); keep only z.txt in the log."]; Abort[]
];
zLines = Take[Select[Drop[lines, First[hdrIdx]], # =!= "" &], UpTo[nExpected]];
If[Length[zLines] < nExpected, Print["Vector is truncated: ", Length[zLines], " of ", nExpected, " entries."]; Abort[]];
z = toNum /@ zLines;

readObj[key_] := Module[{l = Select[lines, StringStartsQ[#, key] &]},
  If[l === {}, Missing["NotInLog"], toNum[Last@StringSplit[First[l], "="]]]];
primalObj = readObj["primalObjective"];
dualObj   = readObj["dualObjective"];

Print["Read ", Length[z], " components from ", logPath,
      ";  precision of z = ", Round[Min[Precision /@ z]], " digits"];


(* ================= Step 2: consistency checks and the bound ================= *)
aF = z[[;; Length[fList]]];
vG = 1/(fList - 1); v2 = 2/(fList + 1); v3 = 1/(fList + 3);

normVal = v3 . aF;
runType = Which[
  Abs[normVal + 1] < 10^-20, "upper",      (* norm = -v3 : g3 <= ... *)
  Abs[normVal - 1] < 10^-20, "lower",      (* norm = +v3 : g3 >= ... *)
  True, Print["v3.a = ", N[normVal, 20], " is neither -1 nor +1: wrong ordering, or not z?"]; Abort[]];
objVal = -(vG + g2Fixed v2) . aF;

Print["v3.a = ", N[normVal, 25], "   (", runType, "-bound run)"];
(* z reproduces the DUAL objective; primal - dual is the duality gap *)
Print["obj.z = ", N[objVal, 25]];
If[! MissingQ[dualObj],   Print["   obj.z - dualObjective   = ", ScientificForm[N[objVal - dualObj, 3]]]];
If[! MissingQ[primalObj], Print["   obj.z - primalObjective = ", ScientificForm[N[objVal - primalObj, 3]]]];
Print["v2.a  = ", N[v2 . aF, 30], "   (slope of the tangent line in the (g2, g3) plane)"];
(* vG.a + g2 v2.a + g3 v3.a >= 0 *)
If[runType === "upper",
  Print["Bound from this functional: g3 <= ", N[vG . aF, 15], " + ", N[v2 . aF, 15], " g2;  at g2 = ",
        g2Fixed, ": g3 <= ", N[vG . aF + g2Fixed v2 . aF, 15]],
  Print["Bound from this functional: g3 >= ", N[-vG . aF, 15], " - ", N[v2 . aF, 15], " g2;  at g2 = ",
        g2Fixed, ": g3 >= ", N[-(vG . aF + g2Fixed v2 . aF), 15]]];
Print["fhat(0) = int f = ", N[Total[aF/(fList + 1)], 15]];


(* =============== Step 3: kernels and cached smearing (= generator) =============== *)
d = 6;
calP[J_, zz_] := Hypergeometric2F1[-J, J + d - 3, (d - 2)/2, (1 - zz)/2];
kappa[J_] := Together[(2 + v)/(1 + v)^2 calP[J, 1 + 2 v] -
    v^2 ((4 + 3 v)/(1 + v)^2 calP[J, 1] + (4 v)/(1 - v^2) (D[calP[J, zz], zz] /. {zz -> 1}))];
chi[k_, J_] := Together[Module[{w},
    (2 + v) calP[J, 1 + 2 v]/(v (1 + v))/(v (1 + v))^(k/2) -
     SeriesCoefficient[(2 + w) (1 - w) (1 + 2 w)/((v - w) w (1 - v) (1 + w) (1 + v + w))*
       calP[J, 1 + 2 w]/(w (1 + w))^(k/2), {w, 0, -1}]]];
QR[expr_] := Module[{num = Numerator[expr], den = Denominator[expr], R},
  If[FreeQ[den, v], {CoefficientList[expr, v], 0},
   R = Together[(1 - v) expr] /. v -> 1;
   {CoefficientList[PolynomialQuotient[num - R den/(1 - v), den, v], v], R}]];

polyMat[{qc_, R_}, nList_] :=
  {Table[qc[[j + 1]] (-1)^j/(n + 2 j + 1), {n, nList}, {j, 0, Length[qc] - 1}], R};
kernelData[J_] := <|
  "F"  -> polyMat[QR[kappa[J]], fList],
  "X4" -> polyMat[QR[chi[4, J]], iList4],
  "X6" -> polyMat[QR[chi[6, J]], iList6]|>;
(* I_n(y) = int_0^1 p^n/(1+y p^2); downward recurrence, stable for y <= 1 *)
Ivals[nList_, y_] := Module[{desc = Reverse[Sort[nList]], tab = <||>},
  Scan[(tab[#] = Hypergeometric2F1[1, (# + 1)/2, (# + 3)/2, -y]/(# + 1)) &, Take[desc, 2]];
  Scan[If[! KeyExistsQ[tab, #], tab[#] = 1/(# + 1) - y tab[# + 2]] &, Drop[desc, 2]];
  tab /@ nList];
evalBlock[{M_, R_}, y_, Iv_] := Module[
  {dd = If[M === {} || First[M] === {}, 0, Last[Dimensions[M]]]},
  If[dd == 0, ConstantArray[0, Length[Iv]], M . (y^Range[0, dd - 1])] + R Iv];
stateRow[kd_, m2_] := Module[{y = 1/m2, IF, IX},
  IF = Ivals[fList, y]; IX = Ivals[iList4, y];
  Join[evalBlock[kd["F"], y, IF]/m2^2,
       evalBlock[kd["X4"], y, IX] m2^(-7),
       evalBlock[kd["X6"], y, Take[IX, Length[iList6]]] m2^(-10)]];

(* m^4-normalised functional; m^4 F -> 2 fhat(b) as m -> infinity at fixed b = 2J/m *)
m2Of[t_] := 1/Sin[t/2]^2;                 (* theta in (0, Pi]; theta = Pi is m = 1 *)
thetaOf[m2_] := 2 ArcSin[1/Sqrt[m2]];
valueAt[kd_, m2_] := m2^2 (stateRow[kd, m2] . z);

golden[f_, a0_, b0_, iters_, prec_] := Module[
  {r = N[(Sqrt[5] - 1)/2, prec], a = a0, b = b0, c, dd, fc, fd},
  c = b - r (b - a); dd = a + r (b - a); fc = f[c]; fd = f[dd];
  Do[If[fc < fd,
        b = dd; dd = c; fd = fc; c = b - r (b - a); fc = f[c],
        a = c; c = dd; fc = fd; dd = a + r (b - a); fd = f[dd]], {iters}];
  If[fc < fd, {c, fc}, {dd, fd}]];

(* refine the local minima of a sampled function; xs sorted, vs = f /@ xs.
   An interior minimum is bracketed by its neighbours; an endpoint minimum
   (e.g. theta = Pi, the threshold m = 1) is bracketed by the endpoint and its
   neighbour, and the endpoint value itself is kept as a candidate. *)
refineMinima[f_, xs_, vs_, below_, iters_, prec_] := Module[{n = Length[xs], idx},
  idx = Select[Range[n], vs[[#]] < below &&
       (# == 1 || vs[[#]] <= vs[[# - 1]]) && (# == n || vs[[#]] <= vs[[# + 1]]) &];
  Map[Function[i, Module[{lo = xs[[Max[i - 1, 1]]], hi = xs[[Min[i + 1, n]]], g},
      g = golden[f, lo, hi, iters, prec];
      First@SortBy[{g, {xs[[i]], vs[[i]]}}, Last]]], idx]];


(* ======================= Step 4: scan of finite states ======================= *)
scanSpin[J_] := Module[{kd = kernelData[J], dt, ths, f, vals, mins, gi, prc},
  dt = dphiScan/(J + 3/2);
  ths = Union[N[Range[dt, Pi, dt], scanPrec], N[{Pi}, scanPrec],
              N[thetaOf[10^Range[logStep, logM2Max, logStep]], scanPrec]];
  f = Function[t, valueAt[kd, m2Of[t]]];
  vals = f /@ ths;
  prc = Min[Accuracy /@ vals];   (* absolute digits: relative precision is meaningless at a pinned zero *)
  gi = First@Ordering[vals, 1];
  mins = refineMinima[f, ths, vals, polishBelow, goldenIters, scanPrec];
  <|"J" -> J,
    "gridMin" -> {m2Of[ths[[gi]]], vals[[gi]]},
    "minima" -> ({m2Of[#[[1]]], #[[2]]} & /@ mins),          (* {m2, m^4 F} *)
    "accuracy" -> prc|>];

If[$KernelCount == 0, LaunchKernels[]];   (* safe to re-evaluate in a notebook *)
DistributeDefinitions[z, fList, iList4, iList6, d, calP, kappa, chi, QR, polyMat, kernelData,
  Ivals, evalBlock, stateRow, m2Of, thetaOf, valueAt, golden, refineMinima, scanSpin,
  dphiScan, logM2Max, logStep, scanPrec, polishBelow, goldenIters];

Print["Scanning ", Length[scanSpins], " spins ..."];
tScan = AbsoluteTiming[
   spinRes = SortBy[ParallelMap[scanSpin, Reverse[scanSpins], Method -> "FinestGrained"], #["J"] &];
][[1]];
Print["  done in ", Round[tScan], " s"];


(* ======================= Step 5: impact-parameter scan ======================= *)
fhat[b_] := Sum[aF[[s]] HypergeometricPFQ[{(fList[[s]] + 1)/2}, {2, (fList[[s]] + 3)/2}, -b^2/4]/(fList[[s]] + 1),
   {s, Length[fList]}];
DistributeDefinitions[aF, fhat];
bs = N[Range[bStep, bScanMax, bStep], bPrec];
fvals = ParallelMap[fhat, bs, Method -> "CoarsestGrained"];
bMins = refineMinima[fhat, bs, fvals, 10^-4, goldenIters, bPrec];
bViol = Select[bMins, #[[2]] < -violTol &];


(* ============================ Step 6: report ============================ *)
statesViol = Flatten[Function[r, {r["J"], #[[1]], #[[2]]} & /@ Select[r["minima"], #[[2]] < -violTol &]] /@ spinRes, 1];
lowPrec = Select[spinRes, #["accuracy"] < 20 &];

Print[""];
Print["--- finite states (values are m^4 x functional; m^4 F -> 2 fhat(b) at large m) ---"];
Print[TableForm[
   ({#["J"],
      ScientificForm[N[Min[Join[{#["gridMin"][[2]]}, #["minima"][[All, 2]]]], 4], 3],
      Module[{best = First@SortBy[Join[{#["gridMin"]}, #["minima"]], Last]},
        Sequence @@ {N[best[[1]], 8], If[#["J"] > 0, N[2 #["J"]/Sqrt[best[[1]]], 6], "-"]}],
      Count[#["minima"], {_, val_} /; val < -violTol],
      Round[#["accuracy"]]} &) /@ spinRes,
   TableHeadings -> {None, {"J", "min", "at m2", "b = 2J/m", "#viol", "abs. digits"}}]];

Print[""];
Print["--- impact parameter: fhat(b), 0 < b <= ", bScanMax, " ---"];
Print["fhat minima below 1e-4: ", {N[#[[1]], 8], ScientificForm[N[#[[2]], 4], 3]} & /@ bMins];

Print[""];
Print["violating states: ", Length[statesViol], ";  violating b: ", Length[bViol]];
If[statesViol =!= {},
  Print["most negative state: J = ", #[[1]], ", m2 = ", N[#[[2]], 10], ", value = ", N[#[[3]], 4]] &@
    First@SortBy[statesViol, Last]];
If[lowPrec =!= {}, Print["WARNING: fewer than 20 reliable absolute digits for J = ", #["J"] & /@ lowPrec,
                         "; increase scanPrec."]];

Export[reportFile,
  Prepend[Flatten[Function[r, {r["J"], N[#[[1]], 20], N[#[[2]], 10], N[r["accuracy"], 4]} & /@
        Join[{r["gridMin"]}, r["minima"]]] /@ spinRes, 1],
   {"J", "m2", "m4F", "absDigits"}], "CSV"];
Print["per-minimum report written to ", reportFile];


(* ================= Step 7: write / merge extraConstraints.m ================= *)
(* exact rationals: the generator applies N[..., prc] to them, so their precision never caps prc *)
toRat[x_] := Rationalize[SetPrecision[x, 40], 10^-30 Abs[x]];
newStates = Union[{#[[1]], toRat[#[[2]]]} & /@ statesViol];
newB = Union[toRat[#[[1]]] & /@ bViol];

{oldStates, oldB} = Block[{extraStates = {}, extraB = {}},
   If[FileExistsQ[extraFile], Get[extraFile]];
   {extraStates, extraB}];
allStates = Union[oldStates, newStates];
allB = Union[oldB, newB];

Export[extraFile,
  "(* generated by checkPositivity.m from " <> FileNameTake[logPath] <> " on " <> DateString[] <> " *)\n" <>
  "(* " <> ToString[Length[newStates]] <> " new states, " <> ToString[Length[newB]] <> " new b points this round *)\n" <>
  "extraStates = " <> ToString[allStates, InputForm] <> ";\n" <>
  "extraB = " <> ToString[allB, InputForm] <> ";\n", "Text"];
Print["extraConstraints.m: ", Length[allStates], " states (", Length[newStates], " new), ",
      Length[allB], " b points (", Length[newB], " new) -> ", extraFile];
If[newStates === {} && newB === {},
  Print["No violation below -", violTol, " on this scan. Before trusting the bound, rerun with a finer dphiScan and larger scanSpins."]];
