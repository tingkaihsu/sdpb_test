result = 1/(2*Pi) * Integrate[ (20*I*x) / ((x^2+25)^2)*Exp[I*m*x]*Exp[-I*x*t],{x,-Infinity, Infinity}, Assumptions->{ m > 0 }];

Print[result//FullSimplify]


(-m + t)/E^(5*Abs[m - t])