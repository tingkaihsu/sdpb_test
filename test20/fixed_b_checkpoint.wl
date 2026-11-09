(* Run in a fresh kernel: wolframscript -file test20/fixed_b_checkpoint.wl
   Or evaluate Get[".../test20/fixed_b_checkpoint.wl"] in Mathematica.
   Section 4 of the CMRS tutorial compares a sequence at fixed b to a
   separately calculated limit. Here the kernels and weights are test20's. *)
Begin["FixedBCheckpoint`"];
ClearAll["FixedBCheckpoint`*"];
dir = DirectoryName[ExpandFileName[$InputFileName]];
check[label_, condition_] := If[TrueQ[condition], Print["PASS  ", label],
  Print["FAIL  ", label]; Abort[]];

(* 1. Load the production kernel definitions without executing the PMP writer. *)
source = Import[FileNameJoin[{dir, "test20.m"}], "Lines"];
Do[
  With[{lines = Select[source, StringStartsQ[#, prefix] &]},
    check["one source definition for " <> prefix, Length[lines] == 1];
    ToExpression[First[lines]]],
  {prefix, {"NAAAA[", "NBBBB[", "g0AAWeight[", "g0BBWeight[",
    "aaNullIndices =", "bbNullIndices ="}}];
Get[FileNameJoin[{dir, "..", "kernels_AB.m"}]];
mAsq = 1/5; massA = Sqrt[mAsq]; wp = 70;
sAt[j_, b_] := (2 j/b)^2;
abNullIndices = {0, 1}; (* The saved z.txt fixture uses the original two mixed rows. *)

(* 2. Strip the positive phase-space weights and multiply each block by s.
   Along M=2J/b, J(J+7)/s -> t=b^2/4.  The normalized angular functions
   tend to F3(mA b) and F4(mA b); the latter multiplies the spin derivative. *)
Get[FileNameJoin[{dir, "impact_AB.wl"}]];
mixedLimit[b_, parity_, indices_: abNullIndices] := Module[{forward, osc},
  {forward, osc} = Transpose[impactABRow[#, b, massA] & /@ indices];
  {-osc/2, forward - parity osc}];
check["b=0: neutral off-diagonal, even and odd charged endpoints",
  mixedLimit[0, 1] === {{-1/2, 0}, {1, 0}} &&
  mixedLimit[0, -1][[2]] === {3, 0}];

(* Include q=k, interior rows, and derivative orders up to the exported kMax. *)
Do[
  With[{indices = Range[0, 44], b = 2, parity = (-1)^j},
    target = mixedLimit[b, parity, indices];
    s = N[sAt[j, b], wp];
    finite = {s NAABB[indices, s, j, massA]/2, s NABAB[indices, s, j, massA]};
    check["all 45 mixed rows at J=" <> ToString[j],
      Max[Abs[Flatten[finite - N[target, wp]]]] < 1/250]],
  {j, {1024, 1025}}];

(* Each row below is a two-component coefficient vector, before contraction
   with {alpha0,alpha1}. The 1/2 belongs only to the neutral off-diagonal. *)
finiteMixed[j_, b_] := With[{s = N[sAt[j, b], wp]},
  {s NAABB[abNullIndices, s, j, massA]/2, s NABAB[abNullIndices, s, j, massA]}];
mixedSequence[j_, b_] := Join[finiteMixed[j, b], {finiteMixed[j + 1, b][[2]]}];
Do[
  target = N[Join[mixedLimit[b, 1], {mixedLimit[b, -1][[2]]}], wp];
  Print["b = ", b, "   limiting rows {neutral 12, charged even, charged odd}: ", N[target, 9]];
  errors = Table[
    residual = Abs[Flatten[mixedSequence[j, b] - target]];
    Print["  even/odd spins ", {j, j + 1}, "   largest absolute row error ", N[Max[residual], 6]];
    Max[residual], {j, {64, 256, 1024}}];
  check["all six mixed components approach their targets at b=" <> ToString[b],
    Greater @@ errors && Last[errors] < 1/50], {b, {2, 10}}];

(* 3. Extract every elastic row at fixed t. After substituting s=J^2/t,
   the ratio of leading coefficients gives the limit when degrees match.
   A smaller numerator degree gives zero; a larger one contradicts the ansatz. *)
leadingRatio[expr_] := Module[{r = Together[expr], n, d, dn, dd},
  n = Numerator[r]; d = Denominator[r];
  dn = Exponent[n, spin]; dd = Exponent[d, spin];
  Which[dn < dd, 0, dn == dd,
    Cancel[Coefficient[n, spin, dn]/Coefficient[d, spin, dd]],
    True, Print["Unexpected growing elastic row: ", expr]; Abort[]]];
elasticRows[s_, j_] := {s NAAAA[aaNullIndices, s, j, massA],
                         s NBBBB[bbNullIndices, s, j, massA]};
limitRows = Map[leadingRatio, elasticRows[spin^2/t, spin], {2}];
check["first elastic null in each sector is -t/2", limitRows[[All, 1]] === {-t/2, -t/2}];
check["all elastic limits are polynomials in t", AllTrue[Flatten[limitRows], PolynomialQ[#, t] &]];

(* 4. Reassemble the saved upper-bound functional in test20 ordering.
   P_A and P_B include the two coupling weights. The entry H12 already
   contains the factor 1/2. Positive semidefiniteness requires both
   diagonals and the determinant; charged positivity requires both parities. *)
data = Import[FileNameJoin[{dir, "z.txt"}], "Lines"];
coeff = ToExpression[StringReplace[#, {"e" -> "*^", "E" -> "*^"}]] & /@ Rest[data];
check["saved functional dimension and normalization",
  First[StringSplit[First[data]]] == ToString[Length[coeff]] &&
  Length[coeff] == 3 + Total[Length /@ {abNullIndices, aaNullIndices, bbNullIndices}] &&
  coeff[[3]] == -1];
{aB, aA} = coeff[[1 ;; 2]];
{alpha, beta, gamma} = TakeList[Drop[coeff, 3], Length /@ {abNullIndices, aaNullIndices, bbNullIndices}];
polys = {2 aA + beta.limitRows[[1]], 2 aB + gamma.limitRows[[2]]};
check["saved diagonal degrees {13,16}", (Exponent[#, t] & /@ polys) === {13, 16}];
Do[
  {pa, pb} = polys /. t -> b^2/4;
  {h12, even} = N[mixedLimit[b, 1].alpha, wp];
  odd = N[mixedLimit[b, -1][[2]].alpha, wp];
  Print["b = ", b, "   {P_A, P_B, det H, charged even, charged odd} = ",
    N[{pa, pb, pa pb - h12^2, even, odd}, 9]], {b, {0, 2, 10}}];

(* 5. Reproduce two off-grid counterexamples independently at finite J.
   Expected pairs {fixed-b limit, finite-spin value} use the saved z.txt.
   A negative diagonal alone suffices to disprove neutral positivity. *)
Do[
  {b, channel, expected} = example;
  j = 50000; s = sAt[j, b];
  finite = If[channel == 1,
    s (aA g0AAWeight[s] + beta.NAAAA[aaNullIndices, s, j, massA]),
    s (aB g0BBWeight[s] + gamma.NBBBB[bbNullIndices, s, j, massA])];
  pair = {polys[[channel]] /. t -> b^2/4, finite};
  Print["b = ", N[b], "   {limit, finite J=50000} = ", N[pair, 10]];
  check["saved negative diagonal at b=" <> ToString[b],
    Max[Abs[pair - expected]] < 10^-9 && Max[pair] < 0],
  {example, {{137/4, 1, {-1223671542, -1224408927}/10^11},
             {617/20, 2, {-2719676045, -2721460394}/10^12}}}];
Print["All fixed-b checkpoints passed."];
End[];
