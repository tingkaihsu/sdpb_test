(* ::Package:: *)

(* ======================================================================== *)
(* dense_check20.wl  (port of dense_check.wl to test20.m)                   *)
(*                                                                          *)
(* Dense off-grid positivity check of an SDPB functional produced by the    *)
(* PMP in test20.m (or src=).  The functional is read from a run log; every *)
(* kernel and block assembler (Fn, Fc, FA, massList) is taken verbatim from *)
(* that file (its PMP writer is never executed), so this tests the          *)
(* sampling, not a retyped kernel.  test20.m has no phase-space weights:    *)
(* phiAA = phiBB = phiAB = 1 here.  Logs with or without the "N 1" header.  *)
(*                                                                          *)
(* Usage (from the repo root or test20/):                                   *)
(*   wolframscript -file test20/dense_check20.wl test20/0_40/1.log          *)
(*       [neutral|charged|all] [I,II,III,IV] [quick]                        *)
(*       [src=test20.m] [bmax3=400] [bmaxlim=5000] [dphi=16] [fresh=1]      *)
(*   Positional args in that order; key=value args anywhere.                *)
(*   Regions: I  J <= 60 lines        II  62 <= J <= 1000 lines              *)
(*            III J = 1e3..1e6 at fixed b <= bmax3, plus the analytic        *)
(*                J->inf limit on b <= bmaxlim (its non-PSD b windows)       *)
(*            IV  M -> inf at fixed J (symbolic)                             *)
(*   Regions I/II are resumable: finished J lines are appended to           *)
(*   dense_check_out/<log>_<sector>/lines_<reg>.m; rerun the same command   *)
(*   after a crash.  fresh=1 discards them.  dphi=8 halves the cost.        *)
(*                                                                          *)
(* OUTPUT: test20/extra_points_<dir>_<log>.m  (e.g. extra_points_0_40_1.m)  *)
(*   extraNeutralPoints = {{J, M^2}, ...}  extraChargedPoints = {...}       *)
(*   extraImpactB = {b, ...}   -- exact rationals, every one a non-PSD      *)
(*   point of the given functional.  Rebuilt after each sector from all     *)
(*   sector results on disk, so neutral and charged can be separate jobs.   *)
(*   The file's header comment carries the 3-line hook for test20.m.        *)
(*                                                                          *)
(* Test functions (all scale free, all invariant under the phase-space      *)
(* congruence diag(Sqrt[phiAA], Sqrt[phiBB])):                              *)
(*   f_i = F_ii / Sum_k |z_k P_k,ii|       relative diagonal                *)
(*   rho = F12^2 / (F11 F22)               rho <= 1  <=>  det >= 0          *)
(*   mu  = lambda_min of {{f1, c}, {c, f2}},  c = F12/Sqrt[S1 S2]           *)
(*         mu < 0  <=>  the block is not PSD.  mu is what gets minimized.   *)
(*   charged (1x1): mu = F / Sum_k |z_k P_k|                                *)
(* Status: mu < -10^-20 VIOLATION;  -10^-20 <= mu < 0 touching (round-off). *)
(* ======================================================================== *)

$HistoryLength = 0;
FileDir = DirectoryName[ExpandFileName[$InputFileName]];
args0 = Rest[$ScriptCommandLine];
kv = Association[Rule @@@ (StringSplit[#, "=", 2] & /@ Select[args0, StringContainsQ[#, "="] &])];
args = Select[args0, !StringContainsQ[#, "="] &];
srcFile = Lookup[kv, "src", "test20.m"];
If[Length[args] < 1,
  Print["usage: wolframscript -file dense_check.wl <log> [neutral|charged|all] [I,II,III,IV] [quick]"];
  Exit[1]];

logPath = With[{a = args[[1]]},
  Which[FileExistsQ[a], ExpandFileName[a],
        FileExistsQ[FileNameJoin[{FileDir, a}]], FileNameJoin[{FileDir, a}],
        True, Print["log not found: ", a]; Exit[1]]];
sectorArg = If[Length[args] >= 2, args[[2]], "neutral"];
sectors   = Switch[sectorArg, "all", {"neutral", "charged"}, _, {sectorArg}];
regions   = If[Length[args] >= 3 && args[[3]] =!= "quick", StringSplit[args[[3]], ","], {"I", "II", "III", "IV"}];
quick     = MemberQ[args, "quick"];
fresh     = Lookup[kv, "fresh", "0"] === "1";   (* fresh=1: ignore resumable lines_*.m *)

wp        = 60;          (* working precision                                 *)
wpHigh    = 120;         (* re-check precision for the worst points           *)
tolViol   = 10^-20;      (* relative violation threshold                      *)
refineMu  = 1/20;        (* refine every local minimum with mu below this     *)
maxRefine = 40;          (* ... at most this many per J line                  *)
goldenIt  = 30;
dPhiDense = Pi/ToExpression[Lookup[kv, "dphi", "16"]];   (* phase step Pi/16 = 4x denser than the PMP; dphi=8 halves the cost *)
bMinGrid  = 1/100;       (* log grid extends to b = 2J/M = bMinGrid           *)
bMaxIII   = ToExpression[Lookup[kv, "bmax3", "400"]];   (* Region III b range (b > bMaxSupp = 100 is unsampled in test20.m) *)
bMaxLim   = ToExpression[Lookup[kv, "bmaxlim", "5000"]]; (* analytic J->inf limit scanned to here *)

say[x__] := (WriteString["stdout", StringJoin[ToString /@ {x}], "\n"]);
stamp[] := DateString[{"Hour", ":", "Minute", ":", "Second"}];
fmt[x_] := Which[x === Infinity, "Infinity", x === Indeterminate, "NaN",
  NumericQ[x] && x == 0, "0", NumericQ[x], ToString[CForm[SetPrecision[x, 14]]], True, ToString[x]];


(* ------------------------------------------------------------------------ *)
(* 1.  The functional                                                       *)
(* ------------------------------------------------------------------------ *)
lines = StringReplace[#, ";" -> ""] & /@ Import[logPath, "Lines"];
hdr = SelectFirst[lines, StringMatchQ[#, "r = " ~~ __ ~~ "mAsq = " ~~ __] &, Missing[]];
If[MissingQ[hdr], say["no 'r = .., mAsq = ..' header in log"]; Exit[1]];
{rRatio, mAsq} = First @ StringCases[hdr,
   "r = " ~~ r : Except[","] .. ~~ "," ~~ Whitespace ~~ "mAsq = " ~~ m__ :> ToExpression /@ {r, m}];
abLine = SelectFirst[lines, StringStartsQ[#, "AB = "] &, Missing[]];
If[MissingQ[abLine], say["no 'AB = {a,b}' line in log"]; Exit[1]];
abRange = ToExpression[StringTrim[StringDrop[abLine, 5]]];
zStart = First @ FirstPosition[lines, s_String /; StringMatchQ[StringTrim[s], DigitCharacter .. ~~ " 1"], {0}];
If[zStart == 0, (* 0_40-style log: no "N 1" header, z starts on line 1 *)
  nz = LengthWhile[lines, StringMatchQ[StringTrim[#], NumberString ~~ ((("e" | "E") ~~ NumberString) | "")] &];
  If[nz == 0, say["no z vector found"]; Exit[1]],
  nz = ToExpression[First @ StringSplit[lines[[zStart]]]]];

zFull = ToExpression[StringReplace[StringTrim[#], {"e" -> "*^", "E" -> "*^"}]] & /@ lines[[zStart + 1 ;; zStart + nz]];
primalLine = SelectFirst[lines, StringStartsQ[#, "primalObjective"] &, Missing[]];
primalObj = If[MissingQ[primalLine], Missing[],
  ToExpression[StringReplace[StringTrim[Last @ StringSplit[primalLine, "="]], {"e" -> "*^"}]]];


(* ------------------------------------------------------------------------ *)
(* 2.  Production definitions from test20v2.m                               *)
(* ------------------------------------------------------------------------ *)
held = ReadList[FileNameJoin[{FileDir, srcFile}], Hold[Expression]];
isPMPdef[h_] := !FreeQ[h, HoldPattern[PMP2SDP[_, _] := _]];
skipQ[h_] := isPMPdef[h] || !FreeQ[h, HoldPattern[Import | LaunchKernels | $ScriptCommandLine | PMP2SDP]] ||
   !FreeQ[h, HoldPattern[FileDir = _]];
With[{r0 = rRatio, m0 = mAsq},
  Scan[If[!skipQ[#], ReleaseHold[#]] &, held];
  If[{rRatio, mAsq} =!= {r0, m0}, say["ABORT: source overwrote r or mAsq"]; Exit[1]]];
mA = Sqrt[mAsq];
abNullIndices = Range @@ abRange;      (* the log, not the current file, decides the AB rows *)
If[DownValues[phiAA] === {}, (* test20.m blocks carry no phase-space weights *)
  phiAA[_] := 1; phiBB[_] := 1; phiAB[_] := 1];

pmpDef = SelectFirst[held, isPMPdef];
npts = First @ Cases[pmpDef, HoldPattern[npts = v_] :> v, Infinity];
prec = wp;
PositiveMatrixWithPrefactor[_, m_] := m;
ReleaseHold /@ Cases[pmpDef,
  s : HoldPattern[(Set | SetDelayed)[(Fn | Fc | FA)[___] | phiSamples | massSamples | massList[__], _]] :> Hold[s],
  Infinity];

expectedN = 3 + Total[Length /@ {abNullIndices, aaNullIndices, bbNullIndices}];
say["log           : ", logPath];
say["r, mAsq       : ", ToString[rRatio, InputForm], ", ", ToString[mAsq, InputForm], "     AB rows ", abNullIndices];
say["N             : ", nz, "  (expected 3 + ", Length[abNullIndices], " + ", Length[aaNullIndices],
    " + ", Length[bbNullIndices], " = ", expectedN, ")"];
If[nz =!= expectedN || zFull[[3]] =!= -1, say["ABORT: dimension or normalization mismatch"]; Exit[1]];
say["bound a1+r a2 : ", N[zFull[[1]] + rRatio zFull[[2]], 20], "   primalObjective ", N[primalObj, 20]];
If[DownValues[Fn] === {} || !NumberQ[N[mgapSq]], say["ABORT: definitions not loaded"]; Exit[1]];

setZ[p_] := (zz = N[zFull, p]; zzAbs = Abs[zz]);
setZ[wp];
{cB, cA} = zz[[1 ;; 2]];
{cAB, cAA, cBB} = TakeList[Drop[zz, 3], Length /@ {abNullIndices, aaNullIndices, bbNullIndices}];
gapOf["neutral"] = mgapSq; gapOf["charged"] = mcSq;


(* ------------------------------------------------------------------------ *)
(* 3.  Point evaluation                                                     *)
(* ------------------------------------------------------------------------ *)
(* sFu: x times the unweighted entries; these have finite fixed-b limits.   *)
evalPt["neutral", j_, x_] := Module[{P = Fn[j, x, 0], F, S1, S2, f1, f2, c, mu, det},
  F = Map[zz . # &, P, {2}];
  S1 = zzAbs . Abs[P[[1, 1]]]; S2 = zzAbs . Abs[P[[2, 2]]];
  f1 = F[[1, 1]]/S1; f2 = F[[2, 2]]/S2; c = F[[1, 2]]/Sqrt[S1 S2];
  mu = (f1 + f2)/2 - Sqrt[((f1 - f2)/2)^2 + c^2];
  det = F[[1, 1]] F[[2, 2]] - F[[1, 2]]^2;
  <|"x" -> x, "J" -> j, "mu" -> mu, "f1" -> f1, "f2" -> f2,
    "rho" -> If[F[[1, 1]] F[[2, 2]] > 0, F[[1, 2]]^2/(F[[1, 1]] F[[2, 2]]), Infinity],
    "omr" -> If[F[[1, 1]] F[[2, 2]] > 0, det/(F[[1, 1]] F[[2, 2]]), -Infinity],
    "F11" -> F[[1, 1]], "F22" -> F[[2, 2]], "F12" -> F[[1, 2]],
    "sFu" -> x {F[[1, 1]]/phiAA[x], F[[1, 2]]/Sqrt[phiAA[x] phiBB[x]], F[[2, 2]]/phiBB[x]},
    "prec" -> Precision[mu]|>];
evalPt["charged", j_, x_] := Module[{P = Fc[j, x, 0][[1, 1]], F, S},
  F = zz . P; S = zzAbs . Abs[P];
  <|"x" -> x, "J" -> j, "mu" -> F/S, "f1" -> F/S, "f2" -> Indeterminate, "rho" -> Indeterminate,
    "omr" -> Indeterminate, "F11" -> F, "F22" -> Indeterminate, "F12" -> Indeterminate,
    "sFu" -> {x F/phiAB[x]}, "prec" -> Precision[F/S]|>];
bOf[p_] := If[p["J"] == 0, 0, 2 p["J"]/Sqrt[p["x"]]];
status[mu_] := Which[mu < -tolViol, "VIOLATION", mu < 0, "touching", True, "ok"];


(* ------------------------------------------------------------------------ *)
(* 4.  Grids and adaptive refinement                                        *)
(* ------------------------------------------------------------------------ *)
logGrid[a_, b_, perDec_] := With[{n = Ceiling[perDec Log10[N[b/a]]]}, Table[a (b/a)^(k/n), {k, 0, n}]];
thresholdGrid[g_] := Join[{g}, Table[g (1 + 10^-k), {k, 1, 12}]];
phaseGrid[sec_, J_, g_, dphi_] := Module[{nn = J + 7/2, K},
  K = Floor[N[nn thetaSupp[sec, g]/dphi, 30]];
  Table[massSupp[sec, k dphi/nn], {k, 1, K}]];
sortedUnique[l_] := DeleteDuplicates[Sort[N[l, wp]], Abs[#1 - #2] < 10^-45 Abs[#1] &];
lineGrid["I", sec_, J_] := With[{g = gapOf[sec]}, sortedUnique @ Join[
    thresholdGrid[g], logGrid[g, 10^10, 32], phaseGrid[sec, J, g, dPhiDense]]];
lineGrid["II", sec_, J_] := With[{g = gapOf[sec]}, sortedUnique @ Join[
    thresholdGrid[g], logGrid[g, Max[10^10, (2 J/bMinGrid)^2], 32], phaseGrid[sec, J, g, dPhiDense]]];

golden[f_, a0_, b0_, it_] := Module[{a = a0, b = b0, g = N[(Sqrt[5] - 1)/2, wp], c, d, fc, fd},
  c = b - g (b - a); d = a + g (b - a); fc = f[c]; fd = f[d];
  Do[If[fc["mu"] < fd["mu"],
      b = d; d = c; fd = fc; c = b - g (b - a); fc = f[c],
      a = c; c = d; fc = fd; d = a + g (b - a); fd = f[d]], {it}];
  If[fc["mu"] < fd["mu"], fc, fd]];

localMinima[mus_] := With[{n = Length[mus]},
  Select[Range[n], (# == 1 || mus[[#]] <= mus[[# - 1]]) && (# == n || mus[[#]] <= mus[[# + 1]]) &]];
(* every violating local minimum is refined (uncapped): each one becomes an extra sample point *)
pickCandidates[mus_] := With[{idx = localMinima[mus]},
  Union[Select[idx, mus[[#]] < 0 &],
    Take[SortBy[Union[Select[idx, mus[[#]] < refineMu &], Take[SortBy[idx, mus[[#]] &], UpTo[3]]], mus[[#]] &],
      UpTo[maxRefine]]]];

(* one point per non-PSD dip: refined minima with mu < -tolViol, merged when they converge to the same x *)
dipPoints[pts_] := Module[{v = SortBy[Select[pts, #["mu"] < -tolViol &], #["x"] &]},
  If[v === {}, {}, MinimalBy[#, #["mu"] &][[1]] & /@
    Split[v, Abs[Log[#2["x"]/#1["x"]]] < 10^-2 &]]];   (* same dip: within 1% in M^2 (dips are >~10% apart) *)

(* profile value for the heat map: log10(1-rho) (neutral) or log10 mu (charged); +1 marks a violation *)
profVal[p_] := If[p["mu"] < -tolViol, 1.,
  N @ Clip[Log10[Max[If[NumericQ[p["omr"]], p["omr"], p["mu"]], 10^-30]], {-30, 0}]];

scanLine[region_, sec_, J_] := Module[{xs, vals, mus, cand, refined, all, srt},
  xs = lineGrid[region, sec, J];
  vals = evalPt[sec, J, #] & /@ xs;
  mus = #["mu"] & /@ vals;
  cand = pickCandidates[mus];
  refined = Table[With[{lo = Log[xs[[Max[i - 1, 1]]]], hi = Log[xs[[Min[i + 1, Length[xs]]]]]},
      golden[evalPt[sec, J, Exp[#]] &, lo, hi, goldenIt]], {i, cand}];
  all = Join[vals, refined];
  srt = SortBy[all, #["mu"] &];
  <|"region" -> region, "J" -> J, "n" -> Length[all], "minima" -> Length[localMinima[mus]],
    "refined" -> Length[cand], "worst" -> First[srt], "top" -> Take[srt, UpTo[5]],
    "dips" -> dipPoints[Join[refined, vals[[Select[localMinima[mus], mus[[#]] < -tolViol &]]]]],
    "nViol" -> Count[all, p_ /; p["mu"] < -tolViol], "nTouch" -> Count[all, p_ /; -tolViol <= p["mu"] < 0],
    "profile" -> ({N[bOf[#]], profVal[#]} & /@ all)|>];


(* ------------------------------------------------------------------------ *)
(* 5.  Fixed-b (eikonal) limit, J -> infinity  [as in fixed_b_checkpoint]    *)
(* ------------------------------------------------------------------------ *)
Get[FileNameJoin[{FileDir, "impact_AB.wl"}]];
mixedLimit[b_, parity_] := Module[{forward, osc},
  {forward, osc} = Transpose[impactABRow[#, b, mA] & /@ abNullIndices];
  {-osc/2, forward - parity osc}];
leadingRatio[expr_] := Module[{r = Together[expr], n, d, dn, dd},
  n = Numerator[r]; d = Denominator[r]; dn = Exponent[n, spin]; dd = Exponent[d, spin];
  Which[dn < dd, 0, dn == dd, Cancel[Coefficient[n, spin, dn]/Coefficient[d, spin, dd]],
    True, $Failed]];
limitRowsCache := limitRowsCache = Map[leadingRatio,
   {s NAAAA[aaNullIndices, s, spin, mA], s NBBBB[bbNullIndices, s, spin, mA]} /. s -> spin^2/tt, {2}];
(* returns the same keys as evalPt for the limiting x * F_unweighted *)
limitPt["neutral", b_] := Module[{t = b^2/4, rA, rB, mix, A, B, H, S1, S2, f1, f2, c},
  {rA, rB} = limitRowsCache /. tt -> t;
  mix = mixedLimit[b, 1][[1]];
  A = 2 cA + cAA . rA; B = 2 cB + cBB . rB; H = cAB . mix;
  S1 = Abs[2 cA] + Abs[cAA] . Abs[rA]; S2 = Abs[2 cB] + Abs[cBB] . Abs[rB];
  f1 = A/S1; f2 = B/S2; c = H/Sqrt[S1 S2];
  <|"b" -> b, "mu" -> (f1 + f2)/2 - Sqrt[((f1 - f2)/2)^2 + c^2], "f1" -> f1, "f2" -> f2,
    "rho" -> If[A B > 0, H^2/(A B), Infinity], "sFu" -> {A, H, B}|>];
limitPt["charged", b_, par_] := Module[{row = mixedLimit[b, par][[2]], F},
  F = cAB . row; <|"b" -> b, "mu" -> F/(Abs[cAB] . Abs[row]), "sFu" -> {F}|>];


(* b windows where a polynomial in t = b^2/4 is negative, t in [0, infinity) *)
negBWindows[p_] := Module[{r, edges, segs},
  r = Sort[Select[Cases[tt /. NSolve[p == 0, tt, Reals, WorkingPrecision -> 50], _?NumericQ], # > 0 &]];
  edges = Join[{0}, r, {Infinity}];
  segs = Select[Partition[edges, 2, 1],
    (p /. tt -> If[#[[2]] === Infinity, #[[1]] + 1, Mean[#]]) < 0 &];
  If[segs === {}, "none", StringRiffle[("[" <> fmt[2 Sqrt[#[[1]]]] <> ", " <>
      If[#[[2]] === Infinity, "inf", fmt[2 Sqrt[#[[2]]]]] <> "]") & /@ segs, " "]]];

(* ------------------------------------------------------------------------ *)
(* 6.  M -> infinity at fixed J  (Region IV, symbolic in 1/x)                *)
(* ------------------------------------------------------------------------ *)
(* nGeg[n,l,w] = 2F1(-n, n+2l; l+1/2; (1-w)/2) for integer n >= 0; truncating *)
(* in (1-w)/2 ~ 1/x is exact to the order kept.  Valid with symbolic J only  *)
(* when every n = J - k >= 0, i.e. J >= 8; smaller J are expanded exactly.   *)
asymOrder = 2;
unweighted["neutral", j_, x_] := With[{P = Fn[j, x, 0]},
  {P[[1, 1]]/phiAA[x], P[[1, 2]]/Sqrt[phiAA[x] phiBB[x]], P[[2, 2]]/phiBB[x]}];
unweighted["charged", j_, x_] := {Fc[j, x, 0][[1, 1]]/phiAB[x]};
asymCoeffs[sec_, j_] := Module[{rows, ser},
  rows = If[IntegerQ[j], unweighted[sec, j, xs],      (* exact Gegenbauer polynomials *)
    Block[{nGeg},
      nGeg[n_, l_, w_] := Sum[Pochhammer[-n, k] Pochhammer[n + 2 l, k]/(Pochhammer[l + 1/2, k] k!) ((1 - w)/2)^k,
        {k, 0, asymOrder + 3}];
      unweighted[sec, j, xs]]];
  ser = Map[Expand[Normal[Series[xs #, {xs, Infinity, asymOrder}]]] &, rows, {2}];
  Table[Map[Coefficient[#, xs, -k] &, ser, {2}], {k, 0, asymOrder}]];   (* {L0,L1,L2}: entry x row *)

polyNonNegReport[p_, var_, from_, parityStep_] := Module[{pe = Expand[p], roots, lead, deg, bad},
  deg = Exponent[pe, var];
  If[deg <= 0, Return[<|"deg" -> 0, "valAtFrom" -> pe, "ok" -> pe >= 0, "bad" -> {}|>]];
  lead = Coefficient[pe, var, deg];
  roots = Select[var /. NSolve[pe == 0, var, Reals, WorkingPrecision -> wp - 10], # > from &];
  bad = Select[Union @ Flatten[Table[from + parityStep Range[Max[0, Floor[(r - from)/parityStep] - 2],
        Floor[(r - from)/parityStep] + 2], {r, roots}]], (pe /. var -> #) < 0 &];
  <|"deg" -> deg, "lead" -> lead, "rootsAbove" -> roots, "valAtFrom" -> (pe /. var -> from),
    "ok" -> lead > 0 && (pe /. var -> from) >= 0 && bad === {}, "bad" -> bad|>];


(* ------------------------------------------------------------------------ *)
(* 7.  Run                                                                  *)
(* ------------------------------------------------------------------------ *)
stem = FileBaseName[DirectoryName[logPath]] <> "_" <> FileBaseName[logPath];
LaunchKernels[];
say["kernels       : ", $KernelCount, "    start ", stamp[]];
ParallelEvaluate[$HistoryLength = 0];

(* merged, exact (rationalized) extra-sampling file for test20.m, rebuilt from whatever sector *)
(* results exist, so neutral and charged can be run as separate jobs                         *)
ratx[x_] := Rationalize[SetPrecision[x, 40], 10^-18 Abs[x]];
writeMergedExtra[] := Module[{file = FileNameJoin[{FileDir, "extra_points_" <> stem <> ".m"}], rd, nP, cP, bN, txt},
  rd[sec_] := With[{f = FileNameJoin[{FileDir, "dense_check_out", stem <> "_" <> sec, "extra_points.m"}]},
    If[FileExistsQ[f], Get[f], <||>]];
  nP = {#[[1]], ratx[#[[2]]]} & /@ Lookup[rd["neutral"], "points", {}];
  cP = {#[[1]], ratx[#[[2]]]} & /@ Lookup[rd["charged"], "points", {}];
  bN = Union[ratx /@ Select[First /@ Join[Lookup[rd["neutral"], "impact", {}], Lookup[rd["charged"], "impact", {}]], # > 0 &]];
  txt = StringJoin[
    "(* Extra sampling points for test20.m -- generated by dense_check20.wl, ", DateString[], "\n",
    "   functional : ", logPath, "   (AB rows ", ToString[abRange], ", r = ", ToString[rRatio, InputForm], ")\n",
    "   sectors    : neutral ", If[rd["neutral"] === <||>, "NOT RUN", ToString[Lookup[rd["neutral"], "regions"]] <>
                     If[TrueQ[rd["neutral"]["quick"]], " (quick)", ""]],
    ";  charged ", If[rd["charged"] === <||>, "NOT RUN", ToString[Lookup[rd["charged"], "regions"]] <>
                     If[TrueQ[rd["charged"]["quick"]], " (quick)", ""]], "\n",
    "   counts     : ", ToString[Length[nP]], " neutral (J, M^2), ", ToString[Length[cP]], " charged (J, M^2), ",
                     ToString[Length[bN]], " impact b\n",
    "   Each point is a place where that functional fails PSD, so each one is a genuine cut.\n\n",
    "   Hook for test20.m:\n",
    "   (1) after the chargedSpins line:\n",
    "       extraNeutralPoints = {}; extraChargedPoints = {}; extraImpactB = {};\n",
    "       Get[FileNameJoin[{FileDir, \"extra_points_", stem, ".m\"}]];\n",
    "   (2) in PMP2SDP, replace the impactNodes line's result by\n",
    "       impactNodes = Join[impactNodes, extraImpactB];\n",
    "   (3) in pols, right after the Fc ParallelTable line:\n",
    "       Flatten[ N[ ParallelTable[ Fn[p[[1]], p[[2]], x], {p, extraNeutralPoints}], prec] ],\n",
    "       Flatten[ N[ ParallelTable[ Fc[p[[1]], p[[2]], x], {p, extraChargedPoints}], prec] ],\n",
    "*)\n\n",
    "extraNeutralPoints = ", ToString[nP, InputForm], ";\n\n",
    "extraChargedPoints = ", ToString[cP, InputForm], ";\n\n",
    "extraImpactB = ", ToString[bN, InputForm], ";\n"];
  Export[file, txt, "Text"];
  say["extra sampling file: ", file, "   (", Length[nP], " neutral, ", Length[cP], " charged, ", Length[bN], " impact b)"]];

runSector[sec_] := Module[{outDir, lightRes, x2Res, sampleRes, linesI = {}, linesII = {}, res3 = <||>,
    res4 = <||>, lineRes, worstAll, recheck, csv, summary = {}, add, plot, Js, cmp, lim},
  outDir = FileNameJoin[{FileDir, "dense_check_out", stem <> "_" <> sec}];
  If[!DirectoryQ[outDir], CreateDirectory[outDir]];
  add[s__] := (say[s]; AppendTo[summary, StringJoin[ToString /@ {s}]]);
  add["==== ", sec, " sector, ", stem, "   (r = ", rRatio, ", AB rows ", abRange, ", N = ", nz, ")"];

  (* 7a. isolated states and the PMP's own sample points: tests the z read-in *)
  If[sec === "neutral",
    lightRes = Module[{F = Map[zz . # &, FA[JA, mA, 0], {2}]},
      {F, Eigenvalues[F]}];
    x2Res = evalPt["neutral", J2, N[m2sq, wp]];
    add["light A+B block   : matrix ", StringRiffle[Map[fmt, lightRes[[1]], {2}], ", "], "  eigenvalues ",
        StringRiffle[fmt /@ lightRes[[2]], ", "], "   ", status[Min[lightRes[[2]]]/Max[Abs[lightRes[[2]]]]]];
    add["X2 block (1,2)    : mu = ", fmt[x2Res["mu"]], "  rho = ", fmt[x2Res["rho"]],
        "  f1 = ", fmt[x2Res["f1"]], "  f2 = ", fmt[x2Res["f2"]], "   ", status[x2Res["mu"]]]];
  DistributeDefinitions[evalPt, scanLine, lineGrid, golden, pickCandidates, dipPoints, localMinima, profVal, bOf,
    zz, zzAbs, wp, tolViol, refineMu, maxRefine, goldenIt, dPhiDense, bMinGrid, gapOf, massList,
    phiSamples, massSamples, npts, prec];
  Js = If[sec === "neutral", Select[neutralSpins, # <= 1000 &], Select[chargedSpins, # <= 1001 &]];
  If[quick, Js = Take[Js, {1, -1, 8}]];
  sampleRes = ParallelMap[Function[j, MinimalBy[evalPt[sec, j, #] & /@ massList[sec, j], #["mu"] &][[1]]],
    Js, Method -> "FinestGrained"];
  With[{w = MinimalBy[sampleRes, #["mu"] &][[1]]},
    add["PMP sample points : ", Length[Js], " spins, worst mu = ", fmt[w["mu"]], " at J = ", w["J"],
      ", M^2 = ", fmt[w["x"]], "   ", status[w["mu"]],
      "   (expect >= -1e-20: confirms the z vector / kernels match the run)"]];

  (* 7b. Regions I and II: every J line, grid + golden refinement of local minima *)
  Do[
    Js = Switch[{reg, sec},
      {"I", "neutral"}, Range[0, 60, 2],
      {"I", "charged"}, Range[0, 60],
      {"II", "neutral"}, Join[Range[62, 400, 2], Range[404, 1000, 4]],
      {"II", "charged"}, Join[Range[61, 400], Flatten[{#, # + 1} & /@ Range[404, 1000, 4]]]];
    If[quick, Js = If[reg === "I", {0, 10, 30, 60}, {100, 202, 300, 1000}] + If[sec === "charged", 1, 0]];
    (* resumable: each finished J line is appended to lines_<reg>.m; a rerun skips those J *)
    Module[{lf = FileNameJoin[{outDir, "lines_" <> reg <> ".m"}], done = {}, todo, chunk},
      If[fresh && FileExistsQ[lf], DeleteFile[lf]];
      If[FileExistsQ[lf], done = Select[ReadList[lf], AssociationQ]];
      done = Select[done, KeyExistsQ[#, "dips"] && MemberQ[Js, #["J"]] &];
      todo = Complement[Js, #["J"] & /@ done];
      say["  region ", reg, ": ", Length[Js], " J lines (", Length[done], " resumed from disk), start ", stamp[]];
      lineRes = done;
      Do[
        chunk = With[{rg = reg}, ParallelMap[KeyDrop[scanLine[rg, sec, #], "profileRaw"] &, c,
          Method -> "FinestGrained"]];
        Scan[PutAppend[#, lf] &, chunk];
        lineRes = Join[lineRes, chunk];
        say["    ", reg, ": ", Length[lineRes], "/", Length[Js], " J lines   ", stamp[],
          "   dips so far: ", Total[Length[#["dips"]] & /@ lineRes]],
        {c, Partition[Reverse[todo], UpTo[2 Max[$KernelCount, 1]]]}]];
    lineRes = SortBy[lineRes, #["J"] &];
    If[reg === "I", linesI = lineRes, linesII = lineRes];
    With[{w = MinimalBy[lineRes, #["worst"]["mu"] &][[1]]["worst"]},
      add["Region ", reg, " (", Length[Js], " J, ", Total[#["n"] & /@ lineRes], " points, ",
        Total[#["refined"] & /@ lineRes], " minima refined): worst mu = ", fmt[w["mu"]],
        " at J = ", w["J"], ", M^2 = ", fmt[w["x"]], ", b = ", fmt[bOf[w]], ", rho = ", fmt[w["rho"]],
        ", f1 = ", fmt[w["f1"]], ", f2 = ", fmt[w["f2"]], "   ", status[w["mu"]],
        "   violating points: ", Total[#["nViol"] & /@ lineRes], " on ",
        Count[lineRes, l_ /; l["nViol"] > 0], " J lines"]];
    say["  region ", reg, " done ", stamp[]],
    {reg, Intersection[regions, {"I", "II"}]}];

  (* 7c. Region III: J = 10^3..10^6 at fixed b, compared with the analytic limit *)
  If[MemberQ[regions, "III"], Module[{J3, bs, pts, byJ, refined, limVals, pars},
    J3 = If[sec === "neutral", {10^3, 10^4, 10^5, 10^6}, Flatten[{#, # + 1} & /@ {10^3, 10^4, 10^5, 10^6}]];
    bs = N[Join[{1/100, 3/100, 1/10, 1/5}, Table[k Pi/(16 mA), {k, 1, Floor[N[bMaxIII 16 mA/Pi]]}]], wp];
    If[quick, J3 = Take[J3, If[sec === "neutral", 2, 4]]; bs = Take[bs, {1, -1, 8}]];
    say["  region III: ", Length[J3], " spins x ", Length[bs], " b values, start ", stamp[]];
    pts = ParallelTable[evalPt[sec, j, (2 j/b)^2], {j, J3}, {b, bs}, Method -> "FinestGrained"];
    (* refine the 3 lowest local minima in b per spin, plus every violating one *)
    refined = With[{J3v = J3, bsv = bs}, ParallelMap[Function[{jc}, With[{j = J3v[[jc[[1]]]], i = jc[[2]]},
        golden[evalPt[sec, j, (2 j/#)^2] &, bsv[[Max[i - 1, 1]]], bsv[[Min[i + 1, Length[bsv]]]],
          If[j >= 10^5, 16, goldenIt]]]],
      Flatten[Table[{jc, i}, {jc, Length[J3]},
        {i, With[{m = localMinima[#["mu"] & /@ pts[[jc]]]},
           Union[Take[SortBy[m, pts[[jc, #]]["mu"] &], UpTo[3]], Select[m, pts[[jc, #]]["mu"] < 0 &]]]}], 1],
      Method -> "FinestGrained"]];
    pars = If[sec === "neutral", {Null}, {1, -1}];
    limVals = Table[If[sec === "neutral", limitPt[sec, b], limitPt[sec, b, par]], {par, pars}, {b, bs}];
    res3 = <|"J" -> J3, "b" -> bs, "pts" -> pts, "refined" -> refined, "limit" -> limVals|>;
    With[{w = MinimalBy[Join[Flatten[pts], refined], #["mu"] &][[1]]},
      add["Region III (J = ", J3, ", b in (0, ", bMaxIII, "]): worst mu = ", fmt[w["mu"]], " at J = ", w["J"],
        ", b = ", fmt[bOf[w]], ", rho = ", fmt[w["rho"]], "   ", status[w["mu"]]]];
    (* convergence to the eikonal limit: compare normalized x*F_unweighted *)
    (* convergence to the eikonal limit, in the scale-free mu, per b window *)
    cmp = Table[With[{par = If[sec === "neutral", 1, (-1)^J3[[jc]]]},
        With[{lim = limVals[[If[sec === "neutral" || par == 1, 1, 2]]], mus = #["mu"] & /@ pts[[jc]]},
          {J3[[jc]], Sequence @@ Table[Max[Abs[Pick[mus - (#["mu"] & /@ lim), # <= bw & /@ bs]]], {bw, {10, 50, 100, bMaxIII}}]}]],
      {jc, Length[J3]}];
    add["  max |mu(J) - mu_limit(b)| over b <= {10, 50, 100, ", bMaxIII, "}: ",
      StringRiffle[("J=" <> ToString[#[[1]]] <> ": " <> StringRiffle[fmt /@ Rest[#], ", "]) & /@ cmp, ";  "]];
    If[Length[J3] >= 2, Do[With[{a = J3[[k]], b2 = J3[[k + 1]]},
        add["  max |mu(J=", a, ") - mu(J=", b2, ")| over the b grid: ",
          fmt[Max[Abs[(#["mu"] & /@ pts[[k]]) - (#["mu"] & /@ pts[[k + 1]])]]]]],
      {k, 1, Length[J3] - 1, If[sec === "neutral", 1, 2]}]];
    With[{l = MinimalBy[Flatten[limVals], #["mu"] &][[1]]},
      add["  analytic J->infinity limit on the b grid: worst mu = ", fmt[l["mu"]], " at b = ", fmt[l["b"]],
        "   ", status[l["mu"]]]];
    If[sec === "neutral", Module[{rA, rB, pA, pB},
      {rA, rB} = limitRowsCache;
      pA = Expand[2 cA + cAA . rA]; pB = Expand[2 cB + cBB . rB];
      add["  limit diagonals A(t), B(t) are polynomials in t = b^2/4 of degree ", Exponent[pA, tt], ", ", Exponent[pB, tt],
        " (leading coefficients ", fmt[Coefficient[pA, tt, Exponent[pA, tt]]], ", ", fmt[Coefficient[pB, tt, Exponent[pB, tt]]], ")"];
      add["  exact J->infinity diagonal-negative b windows:  A < 0 on ", negBWindows[pA], ";  B < 0 on ", negBWindows[pB]]]];
    (* dense analytic limit out to bMaxLim: list every b window where the limit block is not PSD *)
    Module[{bl, parsL, vals, neg, wins},
      bl = N[Join[Range[1/100, 1, 1/100], Table[k Pi/(64 mA), {k, 1, Floor[N[bMaxLim 64 mA/Pi]]}]], wp];
      parsL = If[sec === "neutral", {Null}, {1, -1}];
      Do[
        vals = If[sec === "neutral", limitPt[sec, #], limitPt[sec, #, par]]["mu"] & /@ bl;
        neg = Pick[bl, Thread[vals < -tolViol]];
        wins = If[neg === {}, {}, {First[#], Last[#]} & /@
           Split[neg, Abs[#2 - #1] < 1.01 Pi/(64 mA) + 10^-6 &]];
        res3[{"limitNegWindows", par}] = <|"wins" -> wins,
          "worst" -> MinimalBy[Transpose[{bl, vals}], Last][[1]],
          (* worst b inside each window (golden on the analytic limit) + nodes across wide windows *)
          "bNew" -> Flatten[Table[With[{lo = w[[1]] - Pi/(64 mA), hi = w[[2]] + Pi/(64 mA),
              lf = If[sec === "neutral", limitPt[sec, #] &, limitPt[sec, #, par] &]},
            Module[{bs = golden[lf, Max[lo, 10^-3], hi, goldenIt]["b"]},
              Select[Union[{bs}, Range[w[[1]], w[[2]], Pi/(16 mA)]], lf[#]["mu"] < 0 &]]], {w, wins}]]|>;
        add["  dense J->inf limit, b in (0, ", bMaxLim, "]", If[par === Null, "", " parity " <> ToString[par]],
          ": worst mu = ", fmt[res3[{"limitNegWindows", par}]["worst"][[2]]], " at b = ",
          fmt[res3[{"limitNegWindows", par}]["worst"][[1]]], ";  non-PSD b windows (", Length[wins], "): ",
          If[wins === {}, "none", StringRiffle[("[" <> fmt[#[[1]]] <> ", " <> fmt[#[[2]]] <> "]") & /@ Take[wins, UpTo[40]], " "]]],
        {par, parsL}]];
    say["  region III done ", stamp[]]]];

  (* 7d. Region IV: x -> infinity at fixed J *)
  If[MemberQ[regions, "IV"], Module[{coeffsSym, Lsym, small, smallRes, symRes, entries, pars},
    say["  region IV: start ", stamp[]];
    DistributeDefinitions[asymCoeffs, unweighted, asymOrder];
    pars = If[sec === "neutral", {1}, {1, -1}];
    coeffsSym = asymCoeffs[sec, Jsym];
    Lsym = Map[zz . # &, coeffsSym, {2}];          (* Lsym[[k+1, entry]] *)
    small = If[sec === "neutral", Range[0, 6, 2], Range[0, 7]];
    smallRes = ParallelMap[Function[j, {j, Map[zz . # &, asymCoeffs[sec, j], {2}]}], small];
    entries = If[sec === "neutral", {"11", "12", "22"}, {"c"}];
    Do[With[{L0 = Lsym[[1]] /. (-1)^Jsym -> par},
      add["Region IV leading x*F (x->inf), symbolic J >= 8", If[sec === "charged", " parity " <> ToString[par], ""],
        ": ", Thread[entries -> N[Chop[Expand[L0], 10^-40], 8]]];
      If[sec === "neutral",
        Module[{d = polyNonNegReport[L0[[1]] L0[[3]] - L0[[2]]^2, Jsym, 8, 2]},
          add["  L0_11 ", polyNonNegReport[L0[[1]], Jsym, 8, 2]["ok"], ", L0_22 ", polyNonNegReport[L0[[3]], Jsym, 8, 2]["ok"],
            ", det L0 ", d["ok"], " (", N[d["valAtFrom"], 8], " at J=8; degree ", d["deg"], ")",
            "   [condition gamma0^2 <= 16 a1 a2 when L0 is J independent]"]],
        add["  L0 >= 0 for all J >= 8 of this parity: ",
          polyNonNegReport[L0[[1]] /. Jsym -> 2 jj + If[par == 1, 0, 1], jj, 4, 1]["ok"]]];
      add["  next order L1(J) (coefficient of 1/x): ", Thread[entries -> N[Chop[Expand[Lsym[[2]] /. (-1)^Jsym -> par], 10^-40], 5]]]],
      {par, pars}];
    Do[With[{L0 = r[[2, 1]]},
        add["  J = ", r[[1]], " exact: L0 = ", N[L0, 8], If[sec === "neutral",
          "  det = " <> fmt[L0[[1]] L0[[3]] - L0[[2]]^2], ""]]], {r, smallRes}];
    res4 = <|"Lsym" -> Lsym, "small" -> smallRes|>;
    say["  region IV done ", stamp[]]]];

  (* 7d'. extra sampling points for the PMP.  Every entry is a point where the functional of    *)
  (*      this log is NOT PSD, so adding it is a genuine cut (it changes the SDP).             *)
  (*   points : {J, M^2} from Regions I/II dips, Region III dips at spins the PMP samples,     *)
  (*            and the large PMP spins evaluated at the worst b of each limit window           *)
  (*   impact : b nodes for physicalImpactBlocks (limit windows + Region III dips)              *)
  Module[{pmpSpins = If[sec === "neutral", neutralSpins, chargedSpins], pts = {}, imp = {}, lw, bigJ, cand},
    pts = Join[pts, Flatten[Function[l, {#["J"], #["x"], #["mu"], l["region"]} & /@ l["dips"]] /@ Join[linesI, linesII], 1]];
    If[res3 =!= <||>,
      With[{d3 = Select[Join[Flatten[res3["pts"]], res3["refined"]], #["mu"] < -tolViol &]},
        pts = Join[pts, {#["J"], #["x"], #["mu"], "III"} & /@ Select[d3, MemberQ[pmpSpins, #["J"]] &]];
        imp = Join[imp, {N[bOf[#]], #["mu"], "III"} & /@ d3]];
      lw = Select[Normal[res3], MatchQ[First[#], {"limitNegWindows", _}] &];
      imp = Join[imp, Flatten[Table[{#, Indeterminate, "limit"} & /@ l[[2]]["bNew"], {l, lw}], 1]];
      (* large sampled spins at the new b: added only where the finite-J block is also negative *)
      bigJ = Select[pmpSpins, # >= 400 &];
      cand = Flatten[Outer[{#1, #2} &, bigJ, DeleteDuplicates[First /@ imp]], 1];
      pts = Join[pts, Select[ParallelMap[With[{p = evalPt[sec, #[[1]], N[(2 #[[1]]/#[[2]])^2, wp]]},
          {#[[1]], p["x"], p["mu"], "III-bigJ"}] &, cand], #[[3]] < -tolViol &]]];
    pts = DeleteDuplicatesBy[SortBy[pts, {#[[1]] &, #[[2]] &}], {#[[1]], Round[Log[#[[2]]], 10^-2]} &];
    imp = DeleteDuplicatesBy[SortBy[imp, First], Round[First[#], 10^-4] &];
    add["EXTRA SAMPLING: ", Length[pts], " (J, M^2) points on ", Length[DeleteDuplicates[First /@ pts]],
      " spins (", StringRiffle[KeyValueMap[#1 <> ": " <> ToString[#2] &, Counts[#[[4]] & /@ pts]], ", "], "),  ",
      Length[imp], " impact b nodes"];
    Put[<|"sector" -> sec, "log" -> logPath, "abRange" -> abRange, "points" -> pts, "impact" -> imp,
        "regions" -> regions, "quick" -> quick|>, FileNameJoin[{outDir, "extra_points.m"}]];
    writeMergedExtra[]];

  (* 7e. re-evaluate the 20 worst points at 120 digits *)
  worstAll = Take[SortBy[Join[Flatten[#["top"] & /@ Join[linesI, linesII]],
      If[res3 === <||>, {}, Join[Flatten[res3["pts"]], res3["refined"]]]], #["mu"] &], UpTo[20]];
  setZ[wpHigh]; DistributeDefinitions[zz, zzAbs];
  recheck = With[{wh = wpHigh}, ParallelMap[{#, evalPt[sec, #["J"], SetPrecision[#["x"], wh]]} &, worstAll]];
  setZ[wp]; DistributeDefinitions[zz, zzAbs];
  add["120-digit recheck of the ", Length[recheck], " worst points: max |mu60 - mu120| = ",
    fmt[Max[Abs[#[[1]]["mu"] - #[[2]]["mu"]] & /@ recheck]], ", min precision of mu at 60 digits = ",
    fmt[Min[#[[1]]["prec"] & /@ recheck]], " digits"];
  Do[add["   J = ", r[[2]]["J"], ", M^2 = ", fmt[r[[2]]["x"]], ", b = ", fmt[bOf[r[[2]]]], ", mu120 = ", fmt[r[[2]]["mu"]],
      ", rho = ", fmt[r[[2]]["rho"]], ", f1 = ", fmt[r[[2]]["f1"]], ", f2 = ", fmt[r[[2]]["f2"]], "  ", status[r[[2]]["mu"]]],
    {r, Take[recheck, UpTo[10]]}];

  (* 7f. outputs *)
  csv = Prepend[Table[With[{w = l["worst"]},
      {l["region"], w["J"], fmt[w["x"]], fmt[bOf[w]], fmt[w["f1"]], fmt[w["f2"]], fmt[w["rho"]], fmt[w["omr"]],
       fmt[w["mu"]], fmt[w["F11"]], fmt[w["F22"]], fmt[w["F12"]], l["nViol"], l["nTouch"], l["n"], status[w["mu"]]}],
      {l, Join[linesI, linesII]}],
    {"region", "J", "M2", "b", "F11_rel", "F22_rel", "rho", "one_minus_rho", "mu", "F11", "F22", "F12",
     "n_violating", "n_touching", "n_points", "status"}];
  Export[FileNameJoin[{outDir, "worst_per_J.csv"}], csv, "CSV"];
  If[res3 =!= <||>, Export[FileNameJoin[{outDir, "region3_fixed_b.csv"}],
    Prepend[Flatten[Table[With[{p = res3["pts"][[jc, i]]}, {p["J"], fmt[bOf[p]], fmt[p["mu"]], fmt[p["rho"]],
        fmt[p["f1"]], fmt[p["f2"]], Sequence @@ (fmt /@ p["sFu"])}], {jc, Length[res3["J"]]}, {i, Length[res3["b"]]}], 1],
      {"J", "b", "mu", "rho", "F11_rel", "F22_rel", "sF_unweighted..."}], "CSV"]];
  If[Join[linesI, linesII] =!= {}, Quiet @ Check[
    plot = Module[{data, viol},
      data = Flatten[Table[{Log10[#[[1]]], l["J"], -#[[2]]} & /@ Select[l["profile"], #[[1]] > 0 &],
          {l, Join[linesI, linesII]}], 1];
      viol = Cases[data, {lb_, j_, v_} /; v < 0 :> {lb, j}];
      Show[ListDensityPlot[Select[data, #[[3]] >= 0 &], InterpolationOrder -> 0,
          ColorFunction -> (ColorData["SunsetColors"][1 - #] &), PlotRange -> {All, All, {0, 30}},
          PlotLegends -> BarLegend[Automatic, LegendLabel -> "-log10(1-\[Rho])"],
          FrameLabel -> {"log10 b,  b = 2J/M", "J"}, ImageSize -> 900,
          PlotLabel -> stem <> " " <> sec <> ":  tightness  -log10(1-\[Rho])  (red = violation)"],
        If[viol === {}, Graphics[{}], ListPlot[viol, PlotStyle -> {Red, PointSize[Small]}]]]];
    Export[FileNameJoin[{outDir, "heatmap_tightness.png"}], plot], say["heat map export failed"]]];
  Export[FileNameJoin[{outDir, "summary.txt"}], StringRiffle[summary, "\n"], "Text"];
  Put[<|"linesWorst" -> (#["worst"] & /@ Join[linesI, linesII]), "region3" -> res3, "region4" -> res4,
      "recheck" -> recheck|>, FileNameJoin[{outDir, "results.m"}]];
  say["outputs in ", outDir, "   ", stamp[]];
];

Scan[runSector, sectors];
CloseKernels[];
