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
          {( 1 + x )^2, ( 1 + x )*( 3 - 2*( J + 1 )*J ), -1/36*(J*(1 + J)*(150 + J*(1 + J)*(-43 + 2*J*(1 + J)))), -1/288*((-3 + J)*J*(1 + J)*(4 + J)*(204 + J*(1 + J)*(-32 + J + J^2))), -1/14400*(J*(1 + J)*(246960 + J*(1 + J)*(-67908 + J*(1 + J)*(4916 + J*(1 + J)*(-155 + 2*J*(1 + J)))))), -1/36*((-2 + J)*J*(1 + J)*(3 + J)*(-49 + 2*J*(1 + J)))}
        }}
      |>],
      {J, Jlist}
    ];
    
    norm = {0, 1, 0, 0, 0, 0};
    obj  = {-1, 0, 0, 0, 0, 0};
    
    WritePmpJson[
      jsonFile,
      SDP[obj, norm, pols],
      prec,
      getAnalyticSampleData
    ]
];


test2SDP["pmp.json", 200];
