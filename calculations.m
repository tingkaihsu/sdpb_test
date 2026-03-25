(* Print[334.1 / 326.6]


(* Define the function to solve for F given Rr *)
f[Rr_] := F * Log[2] * ArcCosh[Exp[Log[2]/F]/2] - (Rr - 1)/(Rr + 1)

(* Example: solve for a measured Rr = 1.05 *)
RrMeasured = 1.02;

(* Solve numerically; initial guess F=0.9 works well *)
sol = FindRoot[f[RrMeasured] == 0, {F, 0.9}, WorkingPrecision -> 10]

(* Extract the value *)
Fvalue = F /. sol[[1]]

(* Print nicely *)
Print["For Rr = ", RrMeasured, ", the correction factor F = ", Fvalue]

d = 20 * 10^(-10)
Ra = 685.2
Rb = 220

Print[ Pi*d / Log[2] * (Ra + Rb) / 2 * Fvalue ] *)



q = 1.6 * 10^(-19)
RH = -4.8 * 10^(-11)
rho = 4.1 * 10^(-8)

Print[1/(q*RH) * 10^(-6)]

Print[RH/rho * 10^(4)]