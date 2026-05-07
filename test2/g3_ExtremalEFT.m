(* Finding the lower bound of g3 *)

<<"../SDPB.m";

(* Main definition *)
test2SDP[jsonFile_, prec_:200] := Module[
    {
        Jlist, pols, norm, obj
    },
    
    Jlist = Range[0, 40, 2];  (* J = 0,2,...,40 *)
    
    (* build one PositiveMatrixWithPrefactor per J *)
    pols = Table[
      PositiveMatrixWithPrefactor[<|
        "prefactor"   -> DampedRational[1, {}, 1/E, x],
        "polynomials" -> {{
          {( 1 - x )^3, ( 1 - x )^4*( 3 - 2*( J + 1 )*J ), -1/36*(J*(1 + J)*(150 + J*(1 + J)*(-43 + 2*J*(1 + J))))*(1-x)^6}
        }}
      |>],
      {J, Jlist}
    ];
    
    (* Find upper bound of g3 *)
    norm = {0, -1, 0};
    obj  = {1, 0, 0};
    
    WritePmpJson[
      jsonFile,
      SDP[obj, norm, pols],
      prec,
      getAnalyticSampleData
    ]
];


test2SDP["pmp.json", 200];
