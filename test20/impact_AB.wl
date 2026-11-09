(* D=10, s=M^2. Row n=k(k+1)/2+q has 0<=q<=k in ABAB_check.m.
   Returns {forward, K} before contraction with the mixed-null coefficients.
   s NABAB -> forward - (-1)^J K; s NAABB -> -K. *)
impactABRow[n_Integer?NonNegative, b_, massA_] := Module[{k, q, t, weight},
  k = Floor[(Sqrt[8 n + 1] - 1)/2];
  q = n - k (k + 1)/2;
  t = b^2/4;
  (* Taylor coefficient of order k of the limiting partial wave.
     The k=0 branch also defines the b=0 endpoint without 0^0. *)
  weight = If[k == 0, 1, t^k/(k! Pochhammer[4, k])];
  {If[q == k, 2 weight, 0],
   If[q == 0, weight Hypergeometric0F1[4 + k, -massA^2 t], 0]}
];
