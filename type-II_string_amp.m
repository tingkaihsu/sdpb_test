(* ::Package:: *)

ClearAll["Global`*"];

Amp[s_,t_] := Module[{u = -s-t},
              (s^2+t^2+u^2)^2 * (1/(s*u*t)) *
              (Gamma[1-s] Gamma[1-t] Gamma[1-u]) /
              (Gamma[1+s] Gamma[1+t] Gamma[1+u])];

Print["s-t symmetry: ", Amp[s,t] == Amp[t,s] ];
Print["s-u symmetry: ", Amp[s,t] == Amp[-s-t,t] ];

Print["low-energy Amp[s,t] = ", Series[Amp[s,t], {t,0,3}, {s,0,3}]//Normal//FullSimplify];

Print["low-energy Amp[s,t] = ", Series[Amp[s,t], {s,0,3}, {t,0,3}]//Normal//FullSimplify];



