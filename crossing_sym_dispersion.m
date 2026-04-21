(* kernel *)

ker[k_, s_, p_] := (2s+3p^2)/(s(s+p^2))*((s+p^2)/s^3)^(k/2);

(* M_low *)

Mlow[s_, u_] := 8*Pi*G*(s*t/u + s*u/t + t*u/s) - g*(1/s+1/t+1/u) - g0 + g2*(s^2+t^2+u^2) + g3(s*t*u) + g4(s^2+t^2+u^2)^2/.{t->-s-u};

res = Residue[ker[2, z, p]*Mlow[z, 1/2*z*(-1+Sqrt[z-3p^2]/Sqrt[z+p^2])], {z, 0}];

Print["(k=2) kernel: ", ker[2,z,p]//FullSimplify];

Print["(k=2) sum rule: ", res//FullSimplify];

Print["(k=4) kernel: ", ker[4,z,p]//FullSimplify];

res = Residue[ker[4, z, p]*Mlow[z, 1/2*z*(-1+Sqrt[z-3p^2]/Sqrt[z+p^2])], {z, 0}];
Print["(k=4) sum rule: ", res//FullSimplify];

e1 = (s-4m^2/3) - (-p^2+p^2/(z^3-1)*(z-1)^3);

Print["Massive case, z in terms of s: ", Solve[e1==0, z]//FullSimplify];
(* -1/2*(4*m^2 - 9*p^2 - 3*s + Sqrt[-4*(4*m^2 - 3*s)^2 + (-4*m^2 + 9*p^2 + 3*s)^2])/(4*m^2 - 3*s) *)

t = -p^2+p^2/(z^3-1)*(z-Exp[2Pi*I/3])^3;

-p^2 - (p^2*((-1)^(2/3) + (4*m^2 - 9*p^2 - 3*s + Sqrt[-4*(4*m^2 - 3*s)^2 + (-4*m^2 + 9*p^2 + 3*s)^2])/(2*(4*m^2 - 3*s)))^3)/(-1 - (4*m^2 - 9*p^2 - 3*s + Sqrt[-4*(4*m^2 - 3*s)^2 + (-4*m^2 + 9*p^2 + 3*s)^2])^3/(8*(4*m^2 - 3*s)^3))

Print["t in terms of s and p:", t/.{z->-1/2*(4*m^2 - 9*p^2 - 3*s + Sqrt[-4*(4*m^2 - 3*s)^2 + (-4*m^2 + 9*p^2 + 3*s)^2])/(4*m^2 - 3*s)}//FullSimplify];

e2 = e1/.{m->0};
Print["Massless case, z in terms of s: ", Solve[e2==0, z]//FullSimplify];
Print["t in terms of s and p: ", t/.{z->-1/2*(3*p^2 + s + Sqrt[9*p^4 + 6*p^2*s - 3*s^2])/s}//FullSimplify];
(* -1/2*(3*p^2 + s + Sqrt[9*p^4 + 6*p^2*s - 3*s^2])/s *)

cal = (1+z^3)/(z(1-z^3))*((1-z^3)^2/(-27*p^4*z^3))^(k/2);

Print["cal: ", cal/.{z->-1/2*(3*p^2 + s + Sqrt[9*p^4 + 6*p^2*s - 3*s^2])/s}//FullSimplify];