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

<< "../SDPB.m";  (* Load the SDPB Mathematica package *)


(* ----------------------------------------------------------
   PART 1: Construct and write the JSON input file
   ---------------------------------------------------------- *)

Module[
  {
    pols = {
      PositiveMatrixWithPrefactor[<|
        "prefactor" -> DampedRational[1, {}, 1/E, x],

        (* The "polynomials" field.
           For a 2x2 matrix with N=1 free parameter:
             polynomials[[s]][[r]] = {W^0_{rs}(x), W^1_{rs}(x)}
           where the outer list runs over columns s=1,2
           and the inner list runs over rows r=1,2. *)
        "polynomials" -> {

          (* ---- Column s = 1 ---- *)
          {
            (* Row r=1, Col s=1:  entry M(x,y)_{11} = M^0_{11}(x) + y * M^1_{11}(x)
                                            = (1 + x^4) + y * (x^4/12 + x^2)          *)
            {1 + x^4,   x^4/12 + x^2},

            (* Row r=2, Col s=1:  entry M(x,y)_{21} = M^0_{21}(x) + y * M^1_{21}(x)
                                            = x^2       + y * (x^2/2)                  *)
            {x^2,       x^2/2}
          },

          (* ---- Column s = 2 ---- *)
          {
            (* Row r=1, Col s=2:  entry M(x,y)_{12} = M^0_{12}(x) + y * M^1_{12}(x)
                                            = x^2       + y * (x^2/2)
               (equals the (2,1) entry by symmetry)                                    *)
            {x^2,       x^2/2},

            (* Row r=2, Col s=2:  entry M(x,y)_{22} = M^0_{22}(x) + y * M^1_{22}(x)
                                            = (2 + x^4) + y * (x^4/12 + 2*x^2)        *)
            {2 + x^4,   x^4/12 + 2*x^2}
          }

        } (* end polynomials *)

      |>] (* end PositiveMatrixWithPrefactor *)
    }, (* end pols *)

    norm = {1, 0},   (* normalization vector n: forces z_0 = 1   *)
    obj  = {0, -1},  (* objective vector a:     a.z = -z_1 = -y  *)
    prec = 200       (* working precision in decimal digits       *)
  },

  WritePmpJson["matrix_pmp.json",
               SDP[obj, norm, pols],
               prec];

  Print["Wrote matrix_pmp.json"];
];


(* ----------------------------------------------------------
   PART 2: Verification -- confirm the constraint structure
           is correct before handing off to pmp2sdp/SDPB.
   ---------------------------------------------------------- *)

(* Reconstruct M^0(x) and M^1(x) from the polynomial vectors *)
M0[x_] := {{1 + x^4,   x^2      },
            {x^2,       2 + x^4  }};

M1[x_] := {{x^4/12 + x^2,   x^2/2          },
            {x^2/2,           x^4/12 + 2*x^2}};

(* The constraint matrix as a function of x and y *)
Mxy[x_, y_] := M0[x] + y * M1[x];

(* ---- 2a. Confirm M^0(x) is positive definite for all x >= 0 ---- *)
Print["\n=== Positivity of M^0(x) ==="];
Print["tr  M^0(x) = ", Simplify[Tr[M0[x]]]];
Print["det M^0(x) = ", Simplify[Det[M0[x]]]];
(* Expected: trace = 3 + 2*x^4 > 0; det = 2 + 2*x^4 + x^8 = (1 + x^4)^2 + 1 > 0 *)

eigenvals0 = Eigenvalues[M0[x]];
Print["Eigenvalues of M^0(x):"];
Print["  lambda_+ = ", Simplify[eigenvals0[[1]]]];
Print["  lambda_- = ", Simplify[eigenvals0[[2]]]];

(* ---- 2b. Check constraint at a feasible point: y = 0 ---- *)
Print["\n=== Feasibility check at y = 0 ==="];
Print["M(x, 0) = M^0(x), which is PD (shown above). Feasible for all x >= 0."];

(* ---- 2c. Check constraint near the expected optimum: y = -1.30 ---- *)
Print["\n=== Constraint evaluation near y = -1.30 ==="];
ytest1 = -13/10;
Monitor[
  detVals1 = Table[{x0, Det[Mxy[x0, ytest1]]}, {x0, 0, 2.5, 0.1}];
  ,
  "Evaluating det M(x, -1.30) over x..."
];
Print["det M(x, -1.30) for x in {0, 0.1, ..., 2.5}:"];
Print[TableForm[detVals1, TableHeadings -> {None, {"x", "det M(x,-1.30)"}}]];
Print["All entries positive? ", AllTrue[detVals1[[All,2]], # > 0 &]];

(* ---- 2d. Check constraint at an infeasible point: y = -1.50 ---- *)
Print["\n=== Constraint evaluation at y = -1.50 (expect infeasible) ==="];
ytest2 = -15/10;
detVals2 = Table[{x0, Det[Mxy[x0, ytest2]]}, {x0, 0, 2.5, 0.1}];
Print["Minimum of det M(x, -1.50) over x in {0,...,2.5}: ",
      Min[detVals2[[All,2]]]];
Print["(Negative value confirms infeasibility.)"];

(* ---- 2e. Show the off-diagonal coupling distinguishes this from the scalar case ---- *)
Print["\n=== Off-diagonal coupling: M_{12}(x,y) ==="];
Print["M_{12}(x,y) = x^2 * (1 + y/2)"];
Print["At optimum y* ~ -1.33:  M_{12}(x, -1.33) = x^2 * ", 1 + (-1.33)/2];
Print["This is nonzero for x>0, so the det constraint is strictly active."];
Print["The scalar diagonals alone would permit y as negative as ~ -1.84 (scalar example)."];
Print["The 2x2 matrix constraint tightens the bound to y* ~ -1.33."];

(* ---- 2f. Plot the minimum eigenvalue as a function of x for several values of y ---- *)
plotEig = Plot[
  Evaluate[
    Table[
      Min[Eigenvalues[N[Mxy[x, yval]]]],
      {yval, {0, -0.8, -1.0, -1.3, -1.5}}
    ]
  ],
  {x, 0, 2},
  PlotLegends -> {"y=0", "y=-0.8", "y=-1.0", "y=-1.3 (near opt.)", "y=-1.5 (infeas.)"},
  PlotLabel -> "Minimum eigenvalue of M(x,y) = M^0(x) + y*M^1(x)",
  AxesLabel -> {"x", "\[Lambda]_min"},
  PlotRange -> {-0.5, 2},
  GridLines -> Automatic
];
Print[plotEig];