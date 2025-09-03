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

(* existing code unchanged *)
polynomial[J_] := {
  (1 + x)*(3 - 2*(J + 1)*J),
  (1 + x)^2,
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


Print[matrix]

(* Main definition *)
test2SDP[jsonFile_, prec_:200] := Module[
    {
        pols = {
            PositiveMatrixWithPrefactor[<|
                "prefactor" -> DampedRational[1, {}, 1/E, x],
                "polynomials" -> matrix
            |>]
        },
        norm = {1, 0, 0},
        obj  = {0, 1, 0}
    },
    WritePmpJson[jsonFile, SDP[obj, norm, pols], prec, getAnalyticSampleData]
];


(* A similar computation to the above, except with nontrivial matrix semidefiniteness constraints *)
testSDPMatrix[jsonFile_, prec_:200] := Module[
    {
        pols = {
            PositiveMatrixWithPrefactor[<|
                "prefactor"->DampedRational[1, {}, 1/E, x],
                "polynomials"->
                {{{1 + x^4, 1 + x^4/12 + x^2}, {x^2,     x/5}},
                 {{x^2,     x/5},              {2 + x^4, x^4/3 + 2*x^2}}}
                 |>],
            PositiveMatrixWithPrefactor[<|
                "prefactor"->DampedRational[1, {}, 1/E, x],
                "polynomials"->
                {{{1 + 3x^4/4, 1 + x^4/12 + x^2}, {x^2,     1/2 + x/5}},
                 {{x^2,     1/2 + x/5},        {2 + 3x^4/5, x^4/3 + 2*x^2}}}
                 |>]
        },
        norm = {1, 0},
        obj  = {0, -1}
    },

    WritePmpJson[jsonFile, SDP[obj, norm, pols], prec(*, getAnalyticSampleData*)]
];

test2SDP["pmp.json", 200];
