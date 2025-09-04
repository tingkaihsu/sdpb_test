<<"../SDPB.m";

polynomial[J_] := {
  (1 + x)^2,
  (1 + x)*(3 - 2*(J + 1)*J),
  2*J*(J + 1)*(J*(J + 1) - 8)
};

polynomialsList = Table[polynomial[J], {J, 0, 40, 2}]; (* length 21 *)

extraTriplet = {0, 0, 2};

extendedList = Append[polynomialsList, extraTriplet];

matrix = Table[
  If[i <= j,
    extendedList[[i + 1]],
    extendedList[[j + 1]]
  ],
  {i, 0, 21}, {j, 0, 21}
];

test2SDP[jsonFile_, prec_:200] := Module[
    {
        pols = {
            PositiveMatrixWithPrefactor[<|
                "prefactor" -> DampedRational[1, {}, 1/E, x],
                "polynomials" -> matrix
            |>]
        },
        norm = {0, 1, 0},
        obj  = {-1, 0, 0}
    },
    WritePmpJson[jsonFile, SDP[obj, norm, pols], prec, getAnalyticSampleData]
];

test2SDP["pmp.json", 200];
