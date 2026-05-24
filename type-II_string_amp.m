ClearAll["Global`*"];

Amp[s_,t_] := Module[{u = -s-t},
              (s^2+t^2+u^2)^2 * (-1/(s/4*u/4*t/4)) *
              (Gamma[1-s/4] Gamma[1-t/4] Gamma[1-u/4]) /
              (Gamma[1+s/4] Gamma[1+t/4] Gamma[1+u/4])];

(* result is g2 = 0 *)
Print["g2 = ", SeriesCoefficient[Amp[s,t], {s,0,2},{t,0,0}]//FullSimplify ]

Print["g3 = ", SeriesCoefficient[Amp[s,t], {s,0,3}, {t,0,0}]//FullSimplify ]

