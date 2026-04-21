(* Numerical setup sec5.1 *)
rho1[z_, z0_] := ( Sqrt[4-z0] - Sqrt[4-z] ) / ( Sqrt[4-z0] + Sqrt[4-z] );

rho2[z_, z0_] := ( Sqrt[9-z0] - Sqrt[9-z] ) / ( Sqrt[9-z0] + Sqrt[9-z] );

(* coefficients in ansatz *)

(* Alpha: fully symmetric -- always call with sorted ascending indices *)
al[a_, b_, c_] := alCoeff @@ Sort[{a, b, c}];
 
(* Gamma: fully symmetric -- same canonical form *)
ga[a_, b_, c_] := gaCoeff @@ Sort[{a, b, c}];
 
(* Beta: symmetric only under (a,b,c) <-> (c,b,a) *)
be[a_, b_, c_] := beCoeff[ Min[a, c],  b,  Max[a, c] ];

validTriples[Nmax_Integer] := validTriples[Nmax] =
  Select[
    Flatten[
      Table[
        {a, b, c},
        {a, 0, Nmax},
        {b, 0, Nmax - a},
        {c, 0, Nmax - a - b}
      ],
      2
    ],
    ( #[[1]] * #[[2]] * #[[3]] == 0 ) &
  ];

(* with extension *)
TAAAA[s_, t_, u_, Nmax_Integer] :=
    ap*( 1/(rho1[s, 4/3]-1) + 1/(rho1[t, 4/3]-1) + 1/(rho1[u, 4/3]-1) ) +
    Total[
        Function[{abc},
            al[ abc[[1]], abc[[2]], abc[[3]] ]
            * r1[s, 4/3]^abc[[1]]
            * r1[t, 4/3]^abc[[2]]
            * r1[u, 4/3]^abc[[3]]
        ] /@ validTriples[Nmax]
    ]


TABAB[s_, t_, u_, Nmax_Integer] :=
    - 2/(s - 1)
    - 2/(u - 1)
    + Total[
        Function[{abc},
            be[ abc[[1]], abc[[2]], abc[[3]] ]
            * r2[s, 2/3]^abc[[1]]
            * r1[t, 2/3]^abc[[2]]
            * r2[u, 2/3]^abc[[3]]
        ] /@ validTriples[Nmax]
    ]

TBBBB[s_, t_, u_, Nmax_Integer] :=
    Total[
        Function[{abc},
            ga[ abc[[1]], abc[[2]], abc[[3]] ]
            * r1[s, 0]^abc[[1]]
            * r1[t, 0]^abc[[2]]
            * r1[u, 0]^abc[[3]]
        ] /@ validTriples[Nmax]
    ]

