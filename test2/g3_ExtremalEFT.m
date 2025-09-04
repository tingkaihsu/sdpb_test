(* Examples *)

<<"../SDPB.m";

(* The following is the example in the manual *)
(* Maximize {a,b}.{0,-1} = -b over {a,b} such that {a,b}.{1,0}=a=1 and 

E^(-x)(a(1+x^4) + b(x^4/12 + x^2)) >= 0 for all x>=0

Equivalently,

1+x^4 + b(x^4/12 + x^2) >= 0 for all x >= 0

The prefactor DampedRational[1,{},1/E,x] doesn't affect the answer,
but does affect the choice of sample scalings and bilinear basis.

*)


polynomial[J_] := {
  (1 + x)^2,
  (1 + x)*(3 - 2*(J + 1)*J),
  2*J*(J + 1)*(J*(J + 1) - 8)
};

polynomialsList = Table[polynomial[J], {J, 0, 40, 2}]; (* length 21 *)

(* extra triplet to add *)
extraTriplet = {0, 0, 2};

(* extended list (length 22) *)
extendedList = Append[polynomialsList, extraTriplet];

(* construct 22×22 symmetric matrix from extendedList *)
matrix = Table[
  If[i <= j,
    extendedList[[i + 1]],
    extendedList[[j + 1]]
  ],
  {i, 0, 21}, {j, 0, 21}
];


(* Print[matrix] *)

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
          { ( 1 + x )^2, ( 1 + x )*( 3 - 2*( J + 1 )*J ), 2 * J * ( J + 1 )*( J*( J + 1 ) - 8 ) }
        }}
      |>],
      {J, Jlist}
    ];
    
    norm = {0, 1, 0};
    obj  = {-1, 0, 0};
    
    WritePmpJson[
      jsonFile,
      SDP[obj, norm, pols],
      prec,
      getAnalyticSampleData
    ]
];


(* existing code unchanged
polynomial2[J_] := {
  (1 + x)^2,
  (1 + x)*(3 - 2*(J + 1)*J),
  1/2,
  2*J*(J + 1)*(J*(J + 1) - 8)
};

polynomialsList2 = Table[polynomial2[J], {J, 0, 40, 2}]; (* length 21 *)

(* extra triplet to add *)
extraTriplet2 = {0, 0, 0, 2};

(* extended list (length 22) *)
extendedList2 = Append[polynomialsList2, extraTriplet2];

(* construct 22×22 symmetric matrix from extendedList *)
matrix2 = Table[
  If[i <= j,
    extendedList2[[i + 1]],
    extendedList2[[j + 1]]
  ],
  {i, 0, 21}, {j, 0, 21}
];

(* A similar computation to the above, except with nontrivial matrix semidefiniteness constraints *)
test3SDP[jsonFile_, prec_:200] := Module[
    {
        pols = {
            PositiveMatrixWithPrefactor[<|
                "prefactor"->DampedRational[1, {}, 1/E, x],
                "polynomials"-> matrix2
                 |>]
        },
        norm = {0, 0, 1, 0},
        obj  = {-1, 10.5, 0, 0}
    },
    WritePmpJson[jsonFile, SDP[obj, norm, pols], prec, getAnalyticSampleData]
]; *)

test2SDP["pmp.json", 200];
