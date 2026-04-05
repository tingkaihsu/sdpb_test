(* e1 = Sqrt[m^2 + (pr^2-pi^2+2I*pr*pi)] + Sqrt[mr^2-I*g*mr+(pr^2-pi^2+2I*pr*pi)] - (Sqrt[m^2 + (pr^2-pi^2-2I*pr*pi)] + Sqrt[mr^2+I*g*mr+(pr^2-pi^2-2I*pr*pi)]);


Print[ Solve[e1 == 0, {pr,pi}, Assumptions -> {m>0, mr>0, g>0, pr∈Reals, pi∈Reals}]//FullSimplify ] *)


s = (Sqrt[m^2 + (pr^2-pi^2+2I*pr*pi)] + Sqrt[mr^2-I*g*mr+(pr^2-pi^2+2I*pr*pi)])^2;


Print[ SeriesCoefficient[s/.{pr->pr/e, pi->pi*e}, {e,0,-2}]//FullSimplify ]


t = (Sqrt[m^2 + (pr^2-pi^2+2I*pr*pi)] - Sqrt[m^2 + (pr^2-pi^2-2I*pr*pi)])^2 - (2*I*pi)^2;

Print[ SeriesCoefficient[t/.{pr->pr/e, pi->pi*e}, {e,0,-2}]//FullSimplify ]

Print[ SeriesCoefficient[t/.{pr->pr/e, pi->pi*e}, {e,0,0}]//FullSimplify ]

Print[ Series[t/.{pr->pr/e, pi->pi*e}, {e,0,2}]//FullSimplify ]