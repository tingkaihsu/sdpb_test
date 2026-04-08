(* check_negative_regions_general.m
   GENERAL VERSION for arbitrary dimension of vec{y} (any length of y.txt)
   and arbitrary basis functions.

   CRITICAL CHANGES per your instructions:
   - Checks ONLY intervals BETWEEN consecutive previous sampling points.
     (The SDP solver already guarantees positivity AT the sample points,
      so negativity can ONLY occur in the open intervals between them.)
   - Identifies exactly which pairs of consecutive samples have negative
     regions inside them ("pairs that need to be resampled AGAIN finely").
   - Works for any Length[y] >= 1 (vec{y} can be long).
   - User configures basisFunctions and previousSamplePoints once per run.
   - Dense sampling ONLY inside each inter-sample interval (efficient).
   - Stopping criterion and suggested refinement exactly as in observation.md
     (median x*, N=10 points per negative sub-interval).

   Run after sdpb finishes: mathematica -script check_negative_regions_general.m
*)

SetDirectory[Directory[]];

yFile = "y.txt";
If[!FileExistsQ[yFile],
  Print["ERROR: y.txt not found in current directory!"];
  Quit[]
];

y = ReadList[yFile, Number];
Print["=== LOADED SOLUTION ==="];
Print["vec{y} (length ", Length[y], "): ", y];

(* ====================== USER CONFIGURATION SECTION ======================
   EDIT THESE TWO ITEMS FOR YOUR SPECIFIC PROBLEM.
   Everything else is fully automatic and general. *)

(* Basis functions g_k(x) such that the functional is
      func[x] = Sum[ y[[k]] * basisFunctions[[k]][x] , {k, 1, Length[y]} ]
   (In your PMP, these are exactly the functions whose values were placed
    into the "polynomials" arrays at each sample point.) *)
basisFunctions = {
  Function[x, 1 + x^4],          (* g_1(x) = f1(x) *)
  Function[x, x^4/12 + x^2]      (* g_2(x) = f2(x) *)
  (* Add more lines here for longer vec{y}, e.g.
     Function[x, x^6 + Sin[x]], ... *)
};

(* Previous sampling points used to generate the CURRENT y.txt
   (copy-paste the exact list from your Tests2.m / numeric_pmp.json).
   The script will sort them and check ONLY the intervals between them. *)
previousSamplePoints = {0.06, 0.57, 1.0, 1.25, 1.50};
(* ======================================================================= *)

If[Length[basisFunctions] != Length[y],
  Print["ERROR: Length of basisFunctions (", Length[basisFunctions],
        ") does NOT match Length[y] (", Length[y], ")!"];
  Print["You must define exactly Length[y] basis functions."];
  Quit[]
];

func[x_?NumericQ] := Sum[y[[k]] * basisFunctions[[k]][x], {k, 1, Length[y]}];

(* Sort the previous samples (just in case) *)
samples = Sort[previousSamplePoints];
Print["Checking ", Length[samples]-1, " intervals between the ", Length[samples],
      " previous sampling points: ", samples];

negativeRegions = {};
totalNegLength = 0.0;
nIntervals = Length[samples] - 1;
nPerInterval = 2000;   (* dense interior sampling per interval *)

For[i = 1, i <= nIntervals, i++,
  xa = samples[[i]];
  xb = samples[[i + 1]];
  If[xa >= xb, Continue[]];

  (* ONLY interior points of this specific interval *)
  dx = (xb - xa) / nPerInterval;
  xGrid = Table[xa + k * dx, {k, 1, nPerInterval - 1}];
  funcVals = func /@ xGrid;

  (* Find contiguous negative segments INSIDE this interval *)
  currentStart = None;
  intervalNegRegions = {};
  For[j = 1, j <= Length[xGrid], j++,
    If[funcVals[[j]] < 0,
      If[currentStart === None, currentStart = xGrid[[j]]],
      If[currentStart =!= None,
        AppendTo[intervalNegRegions, {currentStart, xGrid[[j - 1]]}];
        currentStart = None;
      ]
    ]
  ];
  If[currentStart =!= None,
    AppendTo[intervalNegRegions, {currentStart, xb}]
  ];

  If[Length[intervalNegRegions] > 0,
    Print["\nNEGATIVE REGION DETECTED in interval between samples #", i,
          " and #", i+1, " : [", xa, ", ", xb, "]"];
    suggestedForThisInterval = {};

    For[rr = 1, rr <= Length[intervalNegRegions], rr++,
      reg = intervalNegRegions[[rr]];
      xaa = reg[[1]]; xbb = reg[[2]];
      len = xbb - xaa;
      totalNegLength += len;
      xStar = (xaa + xbb)/2;
      s = len / 10;   (* N = 10 as specified in observation.md *)
      newPts = Table[xStar + k * s, {k, -10, 10}];
      AppendTo[suggestedForThisInterval, newPts];

      Print["  Sub-interval [", xaa, ", ", xbb, "]   length = ", len];
      Print["    Median x* = ", xStar, "   spacing s = ", s];
      Print["    Suggested new sample points for this pair: ", newPts];
    ];

    AppendTo[negativeRegions, <|"interval" -> {xa, xb}, "subRegions" -> intervalNegRegions,
                               "suggested" -> Flatten[suggestedForThisInterval]|>];
  ]
];

Print["\n=== FINAL SUMMARY ==="];
Print["Total negative length (sum over all inter-sample intervals) = ", totalNegLength];

If[totalNegLength < 10^-6,
  Print["CONVERGED: Negative regions are negligible (< 1e-6)."],
  Print["REFINE REQUIRED: Add the suggested points (printed above) to your next sampling list."]
];

(* Plot functional over the span of the previous samples (with small margin) *)
If[Length[samples] >= 2,
  xMinPlot = Min[samples] - 0.05*(Max[samples] - Min[samples]);
  xMaxPlot = Max[samples] + 0.05*(Max[samples] - Min[samples]);
  Export["functional_with_solution.png",
    Plot[func[x], {x, xMinPlot, xMaxPlot},
         PlotRange -> All,
         PlotLabel -> "f(x)·y with current solution\n(only inter-sample intervals were checked)",
         AxesLabel -> {"x", "f(x)·y"},
         GridLines -> {samples, None},
         Epilog -> {Red, PointSize[Medium], Point[Transpose[{samples, func /@ samples}]]}]],
    Print["Plot saved as functional_with_solution.png (red points = previous samples)"]
  ]
];

Print["\nTo continue the adaptive loop:"];
Print["1. Copy ALL suggested new points (from the NEGATIVE REGION sections above)"];
Print["2. Add them to the samplePoints list in Tests2.m"];
Print["3. Re-generate numeric_pmp.json"];
Print["4. Delete old sdp/ folder"];
Print["5. Run pmp2sdp + sdpb again"];
Print["6. Re-run this checker."];

Quit[]