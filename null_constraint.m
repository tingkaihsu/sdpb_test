(* Null Constraint from Double Contour *)

Ansatz[s_, u_] := 8*Pi*G*(s*t/u+s*u/t+t*u/s) - l3*(1/s+1/t+1/u) - l4 + g2*(s^2+t^2+u^2) + g3(s*t*u) + g4(s^2+t^2+u^2)^2/.{t->-s-u};

(* Double Contour *)

integrand[s_, u_, n_, l_]:= 1/(s*u)*(Ansatz[s, u]/(s^(n-l)*u^l)-Ansatz[u, s]/(s^l*u^(n-l)));

(* n = 2, l = 1 *)
temp = Residue[integrand[s, u, 2, 1], {s, Infinity}];

res = Residue[temp, {u, 0}];

Print[res//FullSimplify];