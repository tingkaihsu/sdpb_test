(* ::Package:: *)

Print["Import package and momentum-convention setup"]
<< "X`";
momRep = MandelstamRelations[{p1, p2, -p4, -p3, m, s2^(1/2), m, 
     s3^(1/2)} -> {s, t, u}];
(*Eliminate u in terms of s, t*)
Eliu = {u -> m^2 + s2 + s3 + m^2 - (s + t)};


(* ::Section:: *)
(*Trajectories of the two triangle branch points*)


(*
  stri1 and stri2 are the two solutions

    (s2 +/- Sqrt[s2 (s2 - 4 m^2)])/(2 m^2).

  Initially these two triangle branch points are on the triangle
  second sheet.  As s2 is continued from M^2 + i epsilon to the
  unstable pole M^2 - i M Gamma, they enter the triangle first sheet.

  The sheet of the s2 square root and the sheet occupied by a triangle
  branch point are different notions.  The stri functions label the
  two algebraic branch-point solutions, not Riemann sheets.
*)
ClearAll[
  triangleRootUpperRim, triangleRootAfterCut,
  triangleRootOnPath, stri1, stri2,
  stri1OnPath, stri2OnPath
];

(* Physical-cut convention:
   q[x + i 0] = +Sqrt[x (x - 4 m^2)] for x > 4 m^2. *)
triangleRootUpperRim[si_, m_] := I Sqrt[si] Sqrt[4 m^2 - si];

(* Analytic continuation of the same root through the cut. *)
triangleRootAfterCut[si_, m_] := -triangleRootUpperRim[si, m];

(* Algebraic locations of the two triangle branch points, evaluated
   with the upper-rim root convention.  Their triangle-sheet labels
   are assigned separately along the continuation path below. *)
stri1[si_, m_] :=
  (si + triangleRootUpperRim[si, m])/(2 m^2);

stri2[si_, m_] :=
  (si - triangleRootUpperRim[si, m])/(2 m^2);


(* ::Subsection:: *)
(*Numerical path of s2: edit these parameters for the desired pole*)


m0 = 1.;
resonanceM = 2.10;
resonanceWidth = 0.10;

(* Same convention as M^2 - i M Gamma in unsB2. *)
sPole = resonanceM^2 - I resonanceM resonanceWidth;
threshold = 4 m0^2;
iEpsilon = 10^-5 m0^2;

(* Continue the unstable invariant from the upper rim to its pole.
   t = 0: s2 = M^2 + i epsilon.
   t = 1/2: cross the s2 decay cut at M^2.
   t = 1: s2 = M^2 - i M Gamma. *)
s2Start = resonanceM^2 + I iEpsilon;

ClearAll[s2Path, complexCoordinates];

s2Path[t_?NumericQ] := Piecewise[{
  {
    t + I 0,
    0 <= t <= threshold
  },
  {
    resonanceM^2 + I iEpsilon (1 - 2 t),
    0 <= t <= 1/2
  },
  {
    resonanceM^2 - I resonanceM resonanceWidth (2 t - 1),
    1/2 < t <= 1
  }
}];

(* Continue one common square root. Switching its sign at the cut
   keeps both algebraic triangle branch-point trajectories continuous.
   The triangle points themselves go from triangle sheet II to I. *)
triangleRootOnPath[t_?NumericQ] := If[
  t <= 1/2,
  triangleRootUpperRim[s2Path[t], m0],
  triangleRootAfterCut[s2Path[t], m0]
];

stri1OnPath[t_?NumericQ] :=
  (s2Path[t] + triangleRootOnPath[t])/(2 m0^2);

stri2OnPath[t_?NumericQ] :=
  (s2Path[t] - triangleRootOnPath[t])/(2 m0^2);

complexCoordinates[z_?NumericQ] := N[{Re[z], Im[z]}];

(* This must approach zero as iEpsilon -> 0. *)
rootMatchingError = N[
  triangleRootUpperRim[resonanceM^2 + I iEpsilon, m0] -
   triangleRootAfterCut[resonanceM^2 - I iEpsilon, m0]
];
Print["triangle-root sheet-matching error = ", rootMatchingError];

(* Do not sample exactly on the projected cut, where a sheet label is
   ambiguous. *)
sheetGap = 10^-5;
triangleSheetIITimes = Subdivide[0., 1/2 - sheetGap, 350];
triangleSheetITimes = Subdivide[1/2 + sheetGap, 1., 350];

stri1TriangleSheetIIData =
  complexCoordinates[stri1OnPath[#]] & /@ triangleSheetIITimes;
stri1TriangleSheetIData =
  complexCoordinates[stri1OnPath[#]] & /@ triangleSheetITimes;

stri2TriangleSheetIIData =
  complexCoordinates[stri2OnPath[#]] & /@ triangleSheetIITimes;
stri2TriangleSheetIData =
  complexCoordinates[stri2OnPath[#]] & /@ triangleSheetITimes;

cutRoot = Sqrt[resonanceM^2 (resonanceM^2 - threshold)];
stri1CutImage = {
  (resonanceM^2 + cutRoot)/(2 m0^2),
  0.
};
stri2CutImage = {
  (resonanceM^2 - cutRoot)/(2 m0^2),
  0.
};
stri1PoleImage = complexCoordinates[stri1OnPath[1.]];
stri2PoleImage = complexCoordinates[stri2OnPath[1.]];

ClearAll[makeTriangleTrajectoryPlot];

makeTriangleTrajectoryPlot[
   sheetIIData_, sheetIData_, branchName_String, cutImage_,
   poleImage_] :=
 ListLinePlot[
  {sheetIIData, sheetIData},
  PlotStyle -> {
    Directive[Blue, Thick],
    Directive[Red, Thick]
  },
  PlotLegends -> {
    branchName <> ": triangle sheet II",
    branchName <> ": triangle sheet I"
  },
  Frame -> True,
  Axes -> False,
  FrameLabel -> {
    "Re(" <> branchName <> ")",
    "Im(" <> branchName <> ")"
  },
  PlotLabel -> branchName <> " triangle branch-point trajectory",
  PlotRange -> All,
  Epilog -> {
    Black, PointSize[0.014], Point[First[sheetIIData]],
    Text[
      Style["M^2 + i epsilon", 10, Black],
      First[sheetIIData],
      {-1, -1}
    ],

    Darker[Green], PointSize[0.016], Point[cutImage],
    Text[
      Style["cross s2 cut", 10, Darker[Green]],
      cutImage,
      {-1, 1}
    ],

    Purple, PointSize[0.016], Point[poleImage],
    Text[
      Style["image of second-sheet pole", 10, Purple],
      poleImage,
      {-1, 1}
    ]
  },
  ImageSize -> 650
];

stri1TrajectoryPlot = makeTriangleTrajectoryPlot[
  stri1TriangleSheetIIData,
  stri1TriangleSheetIData,
  "stri1",
  stri1CutImage,
  stri1PoleImage
];

stri2TrajectoryPlot = makeTriangleTrajectoryPlot[
  stri2TriangleSheetIIData,
  stri2TriangleSheetIData,
  "stri2",
  stri2CutImage,
  stri2PoleImage
];

(* Display the two trajectories as two separate plots. *)
Print[stri1TrajectoryPlot];
Print[stri2TrajectoryPlot];



temp = LoopIntegrate[(L . L)^2, 
     L, {L + p1, m}, {L + p1 + p2, m}, {L + p1 + p2 + p3, m}] /. 
    momRep /. Eliu;
int = LoopRefine[temp];

Limit[SeriesCoefficient[int, {s, 0, 2}] // KallenExpand, {t -> 0}]

B2[s2_, s3_, m_] := 
  (-DiscB[s2, m, m] + DiscB[s3, m, m])/(2 (s2 - s3)) // 
   DiscExpand;

unsB2[M_, G_, m_] := 
 FullSimplify[
  Limit[B2[s2, s3, m] // DiscExpand, {
    s2 -> M^2 - I G M,
    s3 -> M^2 + I G M
  }]
 ]


Plot[unsB2[x,0.1,1],{x,0.1,6}]

Plot[unsB2[x,0.01,1],{x,0.1,6}]

Plot[Limit[B2[s2,s3,1],{s3->s2}],{s2,0.1,1.9}]



