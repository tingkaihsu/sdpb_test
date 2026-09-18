(* ::Package:: *)

(* 10D ABAB dispersion-relation checks

   B is massless and A has mass m:

     s + t + u = 2 m^2.

   The fixed-t ABAB dispersion relation is kept separate from the
   independently crossed AA -> BB partial-wave expansion.
*)

ClearAll[
  g, gABAB, vldPairs, MABAB, lowEnergyKernel, G0Check,
  phaseABAB10, phaseAABB10, P10, p10Derivative,
  spinMultiplicity10, zABAB10, zAABB10, uCrossed10,
  subtractionKernel10, spectralABAB10, spectralAABB10,
  dispersionKernelABAB10, forwardZeroSubCheck,
  localABAB10, localAABB10, eftCoefficientABAB10,
  eftCoefficientAABB10, crossingNullIRCheck10,
  DblCntABAB, DblCntAABB, ABABNullWeight10,
  AABBNullWeight10, ABABNullWeightClosed10,
  AABBNullWeightClosed10, nullWeightCrossCheck10,
  validDoubleContourLabelQ, validDoubleContourLabels,
  allowedNullSpinQ10, nullBlockWeight10,
  aabbThresholdExponent10, nullConstraintChecks10
];


(* ---------------------------------------------------------------------- *)
(* Low-energy amplitude                                                    *)
(* ---------------------------------------------------------------------- *)

g[a_, b_] := gABAB @@ Sort[{a, b}];

vldPairs[Nmax_Integer?NonNegative] :=
  Flatten[
    Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}],
    1
  ];

MABAB[s_, t_, m_, Nmax_Integer?NonNegative] :=
  Total[
    Function[ab,
      g[ab[[1]], ab[[2]]]
        (u - m^2)^ab[[1]]
        (s - m^2)^ab[[2]]
    ] /@ vldPairs[Nmax]
  ] /. u -> 2 m^2 - s - t;

lowEnergyKernel[
  s_, t_, m_, s1_, s2_, k_, Nmax_Integer?NonNegative
] :=
  MABAB[s, t, m, Nmax]/
    ((s - s1) ((s - s1) (s - s2))^(k/2));

(* This k=0 residue is an algebraic coefficient check.  It does not
   establish convergence of an unsubtracted dispersion relation. *)
G0Check =
  Series[
    -Residue[
      lowEnergyKernel[s, t, m, m^2, m^2 - t, 0, 10],
      {s, Infinity}
    ],
    {t, 0, 2}
  ];


(* ---------------------------------------------------------------------- *)
(* 10D phase space and partial waves                                       *)
(* ---------------------------------------------------------------------- *)

(* Positive convention-dependent overall constants are omitted. *)
phaseABAB10[s_, m_] := (s - m^2)^7/s^4;

(* Geometric phase factor for the separate AA -> BB block. *)
phaseAABB10[s_, m_] := s^(5/4) (s - 4 m^2)^(7/4);

(* In D=10, nu=(D-3)/2=7/2.  The normalization gives P10[J,1]=1. *)
P10[J_, z_] :=
  GegenbauerC[J, 7/2, z]/GegenbauerC[J, 7/2, 1];

(* Dimension of the rank-J symmetric traceless SO(9) representation. *)
spinMultiplicity10[J_Integer?NonNegative] :=
  (2 J + 7) Factorial[J + 6]/(Factorial[7] Factorial[J]);


(* ---------------------------------------------------------------------- *)
(* Channel kinematics                                                      *)
(* ---------------------------------------------------------------------- *)

zABAB10[s_, t_, m_] :=
  1 + 2 s t/(s - m^2)^2;

zAABB10[s_, t_, m_] :=
  (2 t + s - 2 m^2)/Sqrt[s (s - 4 m^2)];

uCrossed10[s_, t_, m_] :=
  2 m^2 - s - t;


(* ---------------------------------------------------------------------- *)
(* Correct fixed-t ABAB dispersive kernel                                  *)
(* ---------------------------------------------------------------------- *)

subtractionKernel10[x_, s1_, s2_, k_] :=
  1/((x - s1) ((x - s1) (x - s2))^(k/2));

spectralABAB10[s_, t_, m_, J_] :=
  phaseABAB10[s, m] P10[J, zABAB10[s, t, m]];

(* Independent coupled-channel building block.  It is not used in the
   fixed-t ABAB dispersion kernel below. *)
spectralAABB10[s_, t_, m_, J_] :=
  phaseAABB10[s, m] P10[J, zAABB10[s, t, m]];

(* Fixed-t structure:

     [K(s') - K(2 m^2 - s' - t)] Im T_ABAB(s',t).

   Both terms multiply the same physical ABAB absorptive partial wave because of the s-u symmetry.
*)
dispersionKernelABAB10[s_, t_, m_, k_, J_] :=
  With[
    {
      s1loc = m^2,
      s2loc = m^2 - t,
      uloc = uCrossed10[s, t, m]
    },
    spectralABAB10[s, t, m, J] (
      subtractionKernel10[s, s1loc, s2loc, k]
        - subtractionKernel10[uloc, s1loc, s2loc, k]
    )
  ];
  
g0ABAB =
  FullSimplify[
    SeriesCoefficient[
      dispersionKernelABAB10[s, t, m, 0, J],
      {t, 0, 0}
    ],
    Assumptions ->
      m > 0 && s > m^2 &&
      Element[{s, m}, Reals] &&
      Element[J, Integers] && J >= 0
  ];
  g0ABAB/.{s->x}


(* ---------------------------------------------------------------------- *)
(* Crossing and double-contour null constraints                            *)
(* ---------------------------------------------------------------------- *)

(* Local coordinates at the crossing-related expansion points:

     ABAB:       sAB = m^2 + sigma,  tAB = tau;
     AA -> BB:   sAA = sigma,        tAA = m^2 + tau.

   Hence localABAB10[sigma,tau] = localAABB10[tau,sigma].
*)

localABAB10[sigma_, tau_, m_, Nmax_Integer?NonNegative] :=
  MABAB[m^2 + sigma, tau, m, Nmax];

localAABB10[sigma_, tau_, m_, Nmax_Integer?NonNegative] :=
  MABAB[m^2 + tau, sigma, m, Nmax];

(* checking crossing symmetry equation *)
eftCoefficientABAB10[
  m_, Nmax_Integer?NonNegative,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  SeriesCoefficient[
    localABAB10[sigma, tau, m, Nmax],
    {sigma, 0, k - q},
    {tau, 0, q}
  ];

eftCoefficientAABB10[
  m_, Nmax_Integer?NonNegative,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  SeriesCoefficient[
    localAABB10[sigma, tau, m, Nmax],
    {sigma, 0, q},
    {tau, 0, k - q}
  ];

crossingNullIRCheck10[
  m_, Nmax_Integer?NonNegative,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  FullSimplify[
    eftCoefficientABAB10[m, Nmax, k, q] -
    eftCoefficientAABB10[m, Nmax, k, q]
  ];
  
crossingNullIRCheck10[m,10,1,0]
crossingNullIRCheck10[m,10,2,0]


(* ---------------------------------------------------------------------- *)
(* Spectral weights after deforming one contour                            *)
(* ---------------------------------------------------------------------- *)

(* x is a physical spectral invariant, never a small contour variable.

   ABAB is expanded in sigma=sAB-m^2 and tau=tAB.  Deforming the
   sigma contour gives the spectral denominator
   (x-m^2)^(k-q+1).

   Crossing maps the AABB expansion point to (sAA,tAA)=(0,m^2).
   Its remaining contour variable is sigma=tAA-m^2, so its partial
   wave must be evaluated at tAA=m^2+sigma.
*)

DblCntABAB[
  x_, tau_, m_, J_,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  spectralABAB10[x, tau, m, J]/
    ((x - m^2)^(k - q + 1) tau^(q + 1));

DblCntAABB[
  x_, sigma_, m_, J_,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  spectralAABB10[x, m^2 + sigma, m, J]/
    (x^(q + 1) sigma^(k - q + 1));

ABABNullWeight10[
  x_, m_, J_,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  Residue[
    DblCntABAB[x, tau, m, J, k, q],
    {tau, 0}
  ];

AABBNullWeight10[
  x_, m_, J_,
  k_Integer?NonNegative, q_Integer?NonNegative
] /; q <= k :=
  Residue[
    DblCntAABB[x, sigma, m, J, k, q],
    {sigma, 0}
  ];


ABABNullWeight10[x,m,J,4,2]//FullSimplify
AABBNullWeight10[x,m,J,4,2]//FullSimplify

ABABNullWeight10[x,m,J,5,2]//FullSimplify
AABBNullWeight10[x,m,J,5,2]//FullSimplify



