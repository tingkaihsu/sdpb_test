(* ============================================================
   Matrix PMP Example for SDPB 3.0.0
   
   Problem:
     maximize  -y   (equivalently, minimize y)
     over      y in R,
     such that M^0(x) + y * M^1(x) >= 0  (positive semidefinite)
               for all x >= 0.
   
   2x2 symmetric polynomial matrices:
   
     M^0(x) = ( 1 + x^4     x^2     )
              ( x^2       2 + x^4   )
   
     M^1(x) = ( x^4/12 + x^2    x^2/2           )
              ( x^2/2           x^4/12 + 2*x^2   )
   
   PMP parameters (using pmp2sdp convention, Eq. 3.1 of manual):
     N = 1   (one free variable y, with z = (z0, z1) and n.z = 1)
     J = 1   (one constraint block)
     m_j = 2 (2x2 positive semidefiniteness constraint)
     normalization n = {1, 0}  => z0 = 1, z1 = y
     objective a    = {0, -1}  => a.z = -z1 = -y
   
   Polynomial vector encoding (see Listing 1 and Listing 2 of manual):
     polynomials[[s]][[r]] = { Q^0_{1,rs}(x), Q^1_{1,rs}(x) }
   
     where W^n_1(x) = M^n(x), so Q^n_{1,rs}(x) = M^n(x)_{rs}.
   
     s=1 (column 1):
       r=1: { 1 + x^4,  x^4/12 + x^2 }    <- entry (1,1)
       r=2: { x^2,      x^2/2         }    <- entry (2,1)
   
     s=2 (column 2):
       r=1: { x^2,      x^2/2         }    <- entry (1,2)  [= entry (2,1) by symmetry]
       r=2: { 2 + x^4,  x^4/12 + 2*x^2 }  <- entry (2,2)
   
   Compare with the scalar example (Listing 4 of manual):
     polynomials = {{{ 1 + x^4,  x^4/12 + x^2 }}}
     which is polynomials[[s=1]][[r=1]] = {Q^0_{1,11}, Q^1_{1,11}}.
   ============================================================ *)

WritePmpJsonNumerical[
  file_,
  SDP[objective_, normalization_, positiveMatricesWithPrefactors_],
  prec_,
  getSampleDataFn_:Function[<||>]
] := exportJson[
  file,
  <|
    "objective"     -> toJsonNumberArray[objective,     prec],
    "normalization" -> toJsonNumberArray[normalization, prec],
    "PositiveMatrixWithPrefactorArray" ->
      Table[toJsonObject[pmp, prec, getSampleDataFn], {pmp, positiveMatricesWithPrefactors}]
  |>
];

<< "../SDPB.m";  (* Load the SDPB Mathematica package *)


M0[x_?NumericQ] := {
  {1 + x^4, x^2},
  {x^2, 2 + x^4}
};

M1[x_?NumericQ] := {
  {x^4/12 + x^2, x^2/2},
  {x^2/2, x^4/12 + 2 x^2}
};

f1List = {
  Function[{x}, M0[x][[1, 1]]],
  Function[{x}, M1[x][[1, 1]]]
};

Print["f1List = ", f1List];

f2List = {
  Function[{x}, M0[x][[2, 1]]],
  Function[{x}, M1[x][[2, 1]]]
};

f3List = {
  Function[{x}, M0[x][[1, 2]]],
  Function[{x}, M1[x][[1, 2]]]
};

f4List = {
  Function[{x}, M0[x][[2, 2]]],
  Function[{x}, M1[x][[2, 2]]]
};

norm = {1, 0};  (* z0 = 1, z1 = y *)
obj = {0, -1};  (* a.z = -z1 = -y *)

testNumericalSDP[spFile_String, jsonFile_String, prec_:650] := Module[
  {rawLines, spLines, samplePoints, sampleScalings, polsRegular},

  (* --- Read and parse sampling_points.txt --- *)
  rawLines = ReadList[spFile, String];
  spLines  = Select[rawLines,
               StringLength[StringTrim[#]] > 0
               && !StringStartsQ[StringTrim[#], "#"] &];

  If[Length[spLines] == 0,
    Print["ERROR: no sample points found in ", spFile]; Quit[1]];

  samplePoints = SetPrecision[ToExpression /@ spLines, prec];
  Print["Read ", Length[samplePoints], " x sample points from ", spFile];
  Print["  x-points    : ", samplePoints];

  (* Scalings = prefactor DampedRational[1,{},1/E,x] evaluated at xi = e^{-xi} *)
  sampleScalings = SetPrecision[Exp[-#] & /@ samplePoints, prec];

  (* --- Regular blocks: one per (xi, Jj) pair.
     Polynomials nesting:  {{ Table[{fk(xi,Jj)}, {k,3}] }}
       {{ ... }}  \[LeftArrow] JSON levels 1 and 2  (column list / row, each size 1)
       Table[...] \[LeftArrow] JSON level 3: polynomial vector, one entry per fList[[k]]
       {fk(xi,Jj)} \[LeftArrow] JSON level 4: degree-0 coefficient list (1 element) --- *)
  polsRegular = Table[
    NumericalPositiveMatrixWithPrefactor[<|
      "prefactor"      -> DampedRational[1, {}, 1/E, x],
      "samplePoints"   -> {samplePoints[[i]]},
      "sampleScalings" -> {sampleScalings[[i]]},
      (* TWO outer braces = JSON levels 1 & 2.
         Table produces level-3 list: {{f1(xi,Jj)}, {f2(xi,Jj)}, {f3(xi,Jj)}}.
         Each {fList[[k]][...]} is the level-4 coefficient list (degree-0: 1 entry). *)
      "polynomials" -> {
        {  (* column s=1 *)
          Table[{SetPrecision[f1List[[k]][samplePoints[[i]]], prec]}, {k, Length[f1List]}],  (* row r=1: entry (1,1) *)
          Table[{SetPrecision[f2List[[k]][samplePoints[[i]]], prec]}, {k, Length[f2List]}]   (* row r=2: entry (2,1) *)
        },
        {  (* column s=2 *)
          Table[{SetPrecision[f3List[[k]][samplePoints[[i]]], prec]}, {k, Length[f3List]}],  (* row r=1: entry (1,2) *)
          Table[{SetPrecision[f4List[[k]][samplePoints[[i]]], prec]}, {k, Length[f4List]}]   (* row r=2: entry (2,2) *)
        }
      }
    |>],
    {i, Length[samplePoints]}
  ];

  (* Flatten polsRegular (2D Table \[RightArrow] flat list) *)
  WritePmpJsonNumerical[
    jsonFile,
    SDP[obj, norm, Flatten[polsRegular] ],
    prec
  ];
  Print["Wrote PMP JSON to ", jsonFile];
];


Module[{myArgs, spFile, jsonFile, prec_:650},

  myArgs = If[Length[$ScriptCommandLine] >= 2, Rest[$ScriptCommandLine], {}];

  If[Length[myArgs] >= 1,
    spFile   = myArgs[[1]];
    jsonFile = If[Length[myArgs] >= 2, myArgs[[2]], "n_pmp.json"];

    Print["=== text9.m ==="];
    Print["  sample_points = ", spFile];
    Print["  output_json   = ", jsonFile];
    Print["  precision     = ", prec];

    testNumericalSDP[spFile, jsonFile, prec];
    Quit[0],

    (* Loaded with << as a library \[LongDash] do nothing. *)
    Null
  ]
];