(* Selecting Basis *)

Clear[basisFunctions];

(* For D = odd and D >= 7, we can use a simple set of basis functions that capture the essential behavior. *)

basisFunctions[D_, maxPow_] := Module[{list, k},
  (* include p^2, p^3, and p^n in 2 ... N_f *)
  list = Table[p^n , {n, 2, N_f}];
  list
];

(* Precompute Analytic Kernel *)