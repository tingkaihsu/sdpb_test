(* ::Package:: *)

(* ---------------------------------------------------------------------- *)
(* Formal 11D -> 10D KK crossing check for the DBI-like amplitude        *)
(* ---------------------------------------------------------------------- *)
(*
   This file performs only the first-stage algebraic checks:

     1. construct the reduced 10D amplitudes AB -> AB and AA -> BB;
     2. verify their crossing relation exactly;
     3. verify the 20 double-contour IR null labels (k,q), 1 <= k <= 5;
     4. display the corresponding low-energy coefficients.

   No 11D partial-wave positivity is assumed here. Pole residues and the
   10D spectral reconstruction with NABAB/NAABB belong to the next stage.
*)

ClearAll[
  s, t, m, sigma, tau,
  parentDBI11, uAB10, kkInvariantsABAB10, kkInvariantsAABB10,
  amplitudeABAB10, amplitudeAABB10,
  localABAB10, localAABB10,
  logGammaOneSeries, gammaRatioSeries,
  parentDBILowEnergySeries, lowEnergyCoefficient,
  validCrossingLabels, exactZeroQ, RunDBIKKCrossingCheck
];


(* ---------------------------------------------------------------------- *)
(* Formal crossing-symmetric parent amplitude                             *)
(* ---------------------------------------------------------------------- *)

parentDBI11[sh_, th_, uh_] :=
  -Gamma[-th] Gamma[-uh]/Gamma[1 + sh]
  -Gamma[-sh] Gamma[-uh]/Gamma[1 + th]
  -Gamma[-sh] Gamma[-th]/Gamma[1 + uh];


(* ---------------------------------------------------------------------- *)
(* 10D kinematics after circle reduction                                  *)
(* ---------------------------------------------------------------------- *)
(*
   B is the massless zero mode and A is a real KK mode of mass m.
   Both ABAB and AABB obey s+t+u=2 m^2 in ten dimensions.
*)

uAB10[s_, t_, m_] := 2 m^2 - s - t;

kkInvariantsABAB10[s_, t_, m_] := {
  s - m^2,
  t,
  uAB10[s, t, m] - m^2
};

kkInvariantsAABB10[s_, t_, m_] := {
  s,
  t - m^2,
  uAB10[s, t, m] - m^2
};

amplitudeABAB10[s_, t_, m_] :=
  Apply[parentDBI11, kkInvariantsABAB10[s, t, m]];

amplitudeAABB10[s_, t_, m_] :=
  Apply[parentDBI11, kkInvariantsAABB10[s, t, m]];


(* Local variables used by the double-contour null constraints:

     ABAB:      s = m^2 + sigma,  t = tau,
     AA -> BB:  s = sigma,        t = m^2 + tau.
*)

localABAB10[sigma_, tau_, m_] :=
  amplitudeABAB10[m^2 + sigma, tau, m];

localAABB10[sigma_, tau_, m_] :=
  amplitudeAABB10[sigma, m^2 + tau, m];


(* ---------------------------------------------------------------------- *)
(* Stable low-energy expansion                                            *)
(* ---------------------------------------------------------------------- *)
(*
   Expanding the three Gamma-function terms separately at sh=th=uh=0
   produces spurious singularities which cancel only in their sum. The
   identity

     Log[Gamma[1+z]] = -EulerGamma z
       + Sum[(-1)^r Zeta[r] z^r/r, {r,2,...}]

   gives the analytic combined expansion without evaluating Gamma at a
   singular point. Two extra orders are retained because each parent term
   contains two explicit inverse Mandelstam variables.
*)

logGammaOneSeries[z_, maxOrder_Integer?NonNegative] :=
  -EulerGamma z
  + Sum[(-1)^r Zeta[r] z^r/r, {r, 2, maxOrder}];

gammaRatioSeries[
  numeratorArgument1_, numeratorArgument2_, denominatorArgument_,
  maxOrder_Integer?NonNegative, expansionParameter_
] :=
  (Normal @ Series[
    Exp[
      logGammaOneSeries[
        expansionParameter numeratorArgument1,
        maxOrder
      ]
      + logGammaOneSeries[
        expansionParameter numeratorArgument2,
        maxOrder
      ]
      - logGammaOneSeries[
        expansionParameter denominatorArgument,
        maxOrder
      ]
    ],
    {expansionParameter, 0, maxOrder}
  ]) /. expansionParameter -> 1;

parentDBILowEnergySeries[
  shArgument_, thArgument_, maxTotalDegree_Integer?NonNegative
] := Module[
  {
    sh0, th0, uh0, eps, ratioOrder,
    ratioS, ratioT, ratioU,
    regulatedAmplitude, regularizedOnShell, truncatedSeries
  },

  ratioOrder = maxTotalDegree + 2;

  ratioS = gammaRatioSeries[-th0, -uh0, sh0, ratioOrder, eps];
  ratioT = gammaRatioSeries[-sh0, -uh0, th0, ratioOrder, eps];
  ratioU = gammaRatioSeries[-sh0, -th0, uh0, ratioOrder, eps];

  regulatedAmplitude =
    -ratioS/(th0 uh0)
    -ratioT/(sh0 uh0)
    -ratioU/(sh0 th0);

  regularizedOnShell = Cancel @ Together[
    regulatedAmplitude /. uh0 -> -sh0 - th0
  ];

  truncatedSeries = Expand @ Normal @ Series[
    regularizedOnShell /. {
      sh0 -> eps sh0,
      th0 -> eps th0
    },
    {eps, 0, maxTotalDegree}
  ] /. eps -> 1;

  truncatedSeries /. {
    sh0 -> shArgument,
    th0 -> thArgument
  }
];

lowEnergyCoefficient[series_, sigmaPower_, tauPower_] :=
  Coefficient[
    Coefficient[series, sigma, sigmaPower],
    tau,
    tauPower
  ];


(* Ordering agrees with NABAB/NAABB in test18.m:

   (1,0),(1,1),(2,0),(2,1),(2,2),...,(5,5).
*)

validCrossingLabels[maxK_Integer?Positive : 5] :=
  Flatten[Table[{k, q}, {k, 1, maxK}, {q, 0, k}], 1];

exactZeroQ[expression_] :=
  TrueQ[expression === 0] || TrueQ[FullSimplify[expression] === 0];


(* ---------------------------------------------------------------------- *)
(* Complete first-stage check                                             *)
(* ---------------------------------------------------------------------- *)

RunDBIKKCrossingCheck[maxK_Integer?Positive : 5] := Module[
  {
    assumptions, labels, maxDegree, parentSeries,
    invariantSumABAB, invariantSumAABB,
    crossingDifference, localCrossingDifference,
    coefficientRows, allCoefficientNullsPass, summary
  },

  assumptions =
    Element[{s, t, m, sigma, tau}, Reals] && m > 0;

  labels = validCrossingLabels[maxK];
  maxDegree = Max[First /@ labels];

  invariantSumABAB = FullSimplify[
    Total[kkInvariantsABAB10[s, t, m]],
    assumptions
  ];

  invariantSumAABB = FullSimplify[
    Total[kkInvariantsAABB10[s, t, m]],
    assumptions
  ];

  crossingDifference = FullSimplify[
    amplitudeABAB10[s, t, m] - amplitudeAABB10[t, s, m],
    assumptions
  ];

  (* This is the precise local identity used by the current (k,q) labels. *)
  localCrossingDifference = FullSimplify[
    localABAB10[sigma, tau, m]
      - localAABB10[tau, sigma, m],
    assumptions
  ];

  parentSeries = parentDBILowEnergySeries[
    sigma,
    tau,
    maxDegree
  ];

  coefficientRows = MapIndexed[
    Function[{label, position},
      Module[
        {k, q, coefficientABAB, coefficientAABB, difference},

        {k, q} = label;

        coefficientABAB = lowEnergyCoefficient[
          parentSeries,
          k - q,
          q
        ];

        (* Crossing interchanges the two local contour powers. *)
        coefficientAABB = lowEnergyCoefficient[
          parentSeries,
          q,
          k - q
        ];

        difference = FullSimplify[
          coefficientABAB - coefficientAABB
        ];

        <|
          "Index" -> First[position] - 1,
          "Label" -> label,
          "ABABCoefficient" -> coefficientABAB,
          "AABBCoefficient" -> coefficientAABB,
          "Difference" -> difference,
          "Passed" -> exactZeroQ[difference]
        |>
      ]
    ],
    labels
  ];

  allCoefficientNullsPass = And @@ Lookup[coefficientRows, "Passed"];

  summary = <|
    "ParentDimension" -> 11,
    "ReducedDimension" -> 10,
    "NumberOfNullLabels" -> Length[labels],
    "ABABInvariantSum" -> invariantSumABAB,
    "AABBInvariantSum" -> invariantSumAABB,
    "InvariantSumsPass" -> And[
      exactZeroQ[invariantSumABAB],
      exactZeroQ[invariantSumAABB]
    ],
    "AmplitudeCrossingDifference" -> crossingDifference,
    "AmplitudeCrossingPass" -> exactZeroQ[crossingDifference],
    "LocalCrossingDifference" -> localCrossingDifference,
    "LocalCrossingPass" -> exactZeroQ[localCrossingDifference],
    "CoefficientNullsPass" -> allCoefficientNullsPass,
    "AllChecksPass" -> And[
      exactZeroQ[invariantSumABAB],
      exactZeroQ[invariantSumAABB],
      exactZeroQ[crossingDifference],
      exactZeroQ[localCrossingDifference],
      allCoefficientNullsPass
    ]
  |>;

  Print["Formal DBI KK crossing check: 11D -> 10D"];
  Print["Null labels: ", Length[labels], " (k = 1, ..., ", maxK, ")"];
  Print[
    Grid[
      Prepend[
        ({
            Lookup[#, "Index"],
            Lookup[#, "Label"],
            Lookup[#, "ABABCoefficient"],
            Lookup[#, "AABBCoefficient"],
            Lookup[#, "Difference"],
            Lookup[#, "Passed"]
          } & /@ coefficientRows),
        {
          "index", "{k,q}", "ABAB coefficient", "AABB coefficient",
          "ABAB - AABB", "pass"
        }
      ],
      Frame -> All,
      Alignment -> Left
    ]
  ];
  Print["Summary: ", summary];

  <|
    "Summary" -> summary,
    "CoefficientChecks" -> coefficientRows,
    "ParentLowEnergySeries" -> parentSeries
  |>
];


DBIKKCrossingCheckResult = RunDBIKKCrossingCheck[5];
