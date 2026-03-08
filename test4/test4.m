(* --- Paste this into Mathematica and run --- *)
<<"../SDPB.m"; (* or the correct path to your SDPB helper file *)

(* ----------------- User parameters (change as needed) ----------------- *)
D = 6;                (* spacetime dimension used in paper examples *)
Mval = 1.;            (* set M = 1 (units) as in paper's Table 3 plotting *)
Jmax = 40;            (* max spin; even spins only; paper used Jmax up to 42 etc. *)
Jlist = Range[0, Jmax, 2];
mMax = 6.0;           (* m_max^2 limit in units of M^2 — the paper uses values like 6 *)
mSamples = 80;        (* number of sample points in m (>= ~50 recommended) *)
pSamples = 80;        (* discretization of p integral over [0,M] *)
fitDegree = 10;       (* polynomial degree used to approximate m-dependence; adjust *)
basisPowerCut = 8;    (* highest power / number of basis functions to include *)
outJson = "pmp_fromA1A2.json";
prec = 200;           (* JSON precision you want in WritePmpJson *)

(* ----------------- Basis functions f_n(p) from Table 1 (flexible) ------------- *)
(* Table 1 suggests using powers and combinations like p^3 - p^2, etc.
   We'll build a small set that reproduces the spirit of the paper and is easy to extend. *)
Clear[basisFunctions];
basisFunctions[D_, maxPow_] := Module[{list, k},
  (* include p^3 - p^2, p^4 - p^2, p^5 - p^2, then fractional powers for even D as needed *)
  list = {p^3 - p^2, p^4 - p^2, p^5 - p^2};
  (* add higher integer powers up to maxPow *)
  list = Join[list, Table[p^n - p^2, {n, 6, maxPow}]];
  (* for larger D include half-integer powers optionally *)
  If[EvenQ[D] && D >= 6,
    list = Join[list, Table[p^(n + 1/2), {n, 3, Floor[(maxPow - 1)/1]}]];
  ];
  list
];

fBasis = basisFunctions[D, basisPowerCut];
nBasis = Length[fBasis];

Print["Using ", nBasis, " basis functions for f(p)."];

(* ----------------- Gegenbauer / P_J implementation ----------------- *)
(* We will use Mathematica's GegenbauerC to represent P_J (up to an overall normalization).
   For our positivity construction the normalization is not important (we keep it consistent). *)
Clear[Pj, PjAt1, PjPrimeAt1];
lambdaVal = (D - 3)/2;
Pj[J_, z_] := GegenbauerC[J, lambdaVal, z];
PjAt1[J_] := Pj[J, 1.0];
PjPrimeAt1[J_] := (D[Pj[J, z], z] /. z -> 1.0);

(* ----------------- Cimproved2 integrand from eq. (3.8) ----------------- *)
(* We follow the structure of eq. (3.8) (u = -p^2).  *)
Clear[Cimproved2p];
Cimproved2p[p_?NumericQ, m_?NumericQ, J_] := Module[
  {u, z, term1, Pval, P1, P1p, denom1, denom2, extra},
  u = -p^2;
  (* argument for Gegenbauer: 1 + 2u/m^2 *)
  z = 1.0 + 2.0 u/m;
  Pval = Pj[J, z];
  P1 = PjAt1[J];
  P1p = PjPrimeAt1[J];
  denom1 = m*(m + u)^2; (* note m stands for m^2 in paper; here we pass m = m^2/M^2 if chosen below *)
  (* main first term (paper uses m^2 everywhere; here m is m^2 variable) *)
  term1 = (2.0 m + u) * Pval / (m*(m + u)^2);
  (* subtractions collected into 'extra' term (paper eq. (3.8) has a closed form subtraction) *)
  (* we implement the simple form used in the paper: subtract terms proportional to u^2/m^6 times combination *)
  (* This is an approximation faithful to the displayed eq. (3.8). *)
  extra = (u^2 / m^3) * ( ((4.0 m + 3.0 u) * P1)/(m + u)^2 + (4.0 u * P1p)/(m^2 - u^2) );
  (* final integrand: Cimproved2[u] ~ term1 - extra *)
  term1 - extra
];

(* Note: in the paper their m variable stands for m^2. We will set mVar = m^2/M^2 below to avoid confusion. *)

(* ----------------- Discretize p and m and compute integrals ----------------- *)
pGrid = Table[(k - 0.5)*Mval/pSamples, {k, 1, pSamples}]; (* midpoints for simple Riemann integration *)
mGrid = Table[Mval + (k - 0.5)*(mMax - Mval)/mSamples, {k, 1, mSamples}]; (* sample m (we treat m as m^2) *)

(* Precompute numeric Gegenbauer derivatives at 1 for all J to speed up *)
Do[
  PjAt1Cache[J] = N[PjAt1[J], 40];
  PjPrimeAt1Cache[J] = N[PjPrimeAt1[J], 40];
, {J, Jlist}];

(* numeric integrator (simple Riemann sum) for stability in Mathematica; replace with NIntegrate if desired *)
clearAndComputeIntegrals[] := Module[
  {Ivals, weights, pStep, basisVals, m, J, row},
  pStep = Mval/pSamples;
  weights = Table[pStep, {pSamples}];
  Ivals = Association[]; (* keys: {J, idxBasis} -> function of mGrid as list *)
  For[Jindex = 1, Jindex <= Length[Jlist], Jindex++,
    J = Jlist[[Jindex]];
    Print["Computing integrals for J=", J];
    (* for each basis function, compute integral over p of f_n(p) * Cimproved2p(p,m,J) as function of m *)
    For[b = 1, b <= nBasis, b++,
      basisVals = Table[
        (fBasis[[b]] /. p -> pGrid[[ip]]) * Cimproved2p[pGrid[[ip]], mGrid[[jm]], J] ,
        {jm, 1, mSamples}, {ip, 1, pSamples}
      ];
      (* basisVals is mSamples x pSamples matrix -> integrate over p by summing along second index *)
      integralByM = Table[
        Sum[basisVals[[jm, ip]]*weights[[ip]], {ip, 1, pSamples}],
        {jm, 1, mSamples}
      ];
      Ivals[{J, b}] = N[integralByM, 20];
    ];
  ];
  Ivals
];

Print["Starting numeric integral computation (may take a minute) ..."];
Ivals = clearAndComputeIntegrals[];
Print["Integrals computed."];

(* ----------------- For each J and each basis coefficient, fit dependence on x = m^2/M^2 by polynomial -------------- *)
xGrid = (mGrid/Mval)^2; (* x = (m/M)^2 *)

(* build polynomial approximations: for each J produce vector polynomials p_n(x) such that:
   integral(m) ≈ Sum_n coeff_n * p_n(x). We'll represent each basis coefficient's m-dependence by a polynomial in x. *)
polyCoeffsByJ = Association[]; (* key: J -> list of polynomial functions for each basis index b *)

Do[
  Print["Fitting polynomials for J=", J];
  polysForThisJ = Table[
    Module[{vals = Ivals[{J, b}], fit},
      (* least-squares fit polynomial in x of degree fitDegree *)
      fit = Fit[Transpose[{xGrid, vals}], Table[x^k, {k, 0, fitDegree}], x];
      (* convert Fit result to exact polynomial expression *)
      N[fit, 20]
    ],
    {b, 1, nBasis}
  ];
  polyCoeffsByJ[J] = polysForThisJ;
, {J, Jlist}];

Print["Polynomial fits ready (degree = ", fitDegree, ")."];

(* ----------------- Build 'pols' in the same shape as your example PositiveMatrixWithPrefactor --- *)
(* We'll encode a positivity block per J. Each block will be a 1x1 positive matrix (scalar >= 0),
   with the polynomial in x representing the integrand for the linear functional associated with each basis coefficient.
   Since the original SDP format used polynomials in x, this matches that style.
*)
polsForPMP = Table[
  Module[{polynomialsList},
    (* For each basis coefficient create a polynomial entry. The "polynomials" field is a list-of-lists-of-polynomials
       matching the example shape. Here we create a single-row matrix whose entries are the linear combination
       of basis polynomials multiplied by symbolic coefficients a_1, a_2, ... .
       To be concrete we expose those coefficients as "a_b" (the decision variables) inside the polynomial,
       but WritePmpJson expects numeric coefficients — different PMP encodings vary. Instead, we will
       embed the basis polynomials as separate polynomial entries; later the SDP/objective can combine them.
       For compatibility with the example, we'll put the polynomials in a single row (1x(nBasis)) so the JSON
       will contain them and your downstream code building the linear combinations can use them.
    *)
    polynomialsList = { Table[ polyCoeffsByJ[J][[b]] , {b, 1, nBasis} ] };
    PositiveMatrixWithPrefactor[ <|
      "prefactor" -> DampedRational[1, {}, 1/E, x],  (* same pattern as example *)
      "polynomials" -> polynomialsList
    |> ]
  ],
  {J, Jlist}
];

(* ----------------- Objective and normalization ---------------------------------------------------------- *)
(* For constructing the SDP that yields inequalities like the ones in Table 3, one usually sets objective & norm.
   We'll create a dummy objective that aims to, e.g., minimize g2 subject to positivity constraints.
   The exact objective vector depends on your normalization of variables. Here we place placeholders. *)

(* Example objective and norm placeholders:
   - norm picks which low-energy coefficients are used to fix scale (like norm = {0,0,1,0} in your initial example)
   - obj encodes the linear objective in the low-energy coupling vector
   You should adapt obj/norm to match exactly how your WritePmpJson/SDP expects mapping to [8πG, g2, g3, ...]. *)
normVec = {0, 0, 1, 0}; (* example: normalize g2 = 1 *)
objVec = {-1, 0, 0, 0}; (* example: minimize Newton's constant contribution; change to your desired objective *)

(* ----------------- Write PMP JSON file --------------------------------------------------------------- *)
sdpProblem = SDP[objVec, normVec, polsForPMP];

Print["Writing JSON to ", outJson, " ..."];
WritePmpJson[outJson, sdpProblem, prec, getAnalyticSampleData];

Print["Done. The file ", outJson, " was written. You can now feed it to SDPB or inspect the JSON."];
