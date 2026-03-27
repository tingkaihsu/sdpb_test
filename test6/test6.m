ClearAll["Global`*"]

y[m2_, J_, p_] := m2 - p^2 + J (J + 1)/(m2 + p^2)
f[n_, p_] := p^n - p^2

(* Numeric integrator – requires exact or high‑precision inputs *)
intNum[m_?NumericQ, n_?IntegerQ, J_?IntegerQ] := 
 NIntegrate[f[n, p] y[m2, J, p], {p, 0, 1}, 
  WorkingPrecision -> 200, PrecisionGoal -> Automatic, AccuracyGoal -> 200]

(* Exact rational argument *)
Print["intNum[(3/2)^2, 3, 2] = ", intNum[(3/2)^2, 3, 2] // N[#, 2] &]

m2x[x_] := 1/(1-x);

Jmax = 40;




test2SDP[jsonFile_, prec_:200] := Module[
    {
        Jlist, pols, norm, obj
    },
    
    Jlist = Range[0, 40, 2];  (* J = 0,2,...,40 *)
    
    (* build one PositiveMatrixWithPrefactor per J *)
    pols = Table[
      PositiveMatrixWithPrefactor[<|
        (* what is damped rational *)
        (* Only damped-rational prefactors are supported. Such prefactors are relevant to the conformal bootstrap problems. *)
        "prefactor"   -> DampedRational[1, {}, 1/E, x],
        "polynomials" -> {{
          { ( 1 + x )^2, ( 1 + x )*( 3 - 2*( J + 1 )*J ), 1/2, 2 * J * ( J + 1 )*( J*( J + 1 ) - 8 ) }
        }}
      |>],
      {J, Jlist}
    ];
    
    norm = {0, 0, 1, 0};
    obj  = {-1, 10.5, 0, 0};
    
    WritePmpJson[
      jsonFile,
      SDP[obj, norm, pols],
      prec,
      getAnalyticSampleData
    ]
];

