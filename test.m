Print["st coefficient = ", SeriesCoefficient[SeriesCoefficient[g2*(s^2+t^2+u^2)+g3*(s*t*u)/.{u -> -s-t}, {s,0,1}], {t,0,1}]//FullSimplify ]
Print["s^2 coefficient = ", SeriesCoefficient[SeriesCoefficient[g2*(s^2+t^2+u^2)+g3*(s*t*u)/.{u -> -s-t}, {s,0,2}], {t,0,0}]//FullSimplify ]


Print["s^3 + t^3 + u^3 = ", s^3+t^3+u^3/.{t -> -s-u}//FullSimplify]

Print["u^3 coeff = ", SeriesCoefficient[-3*s*u*(s + u), {u,0,3}]//FullSimplify]

Print["su^2 coeff = ", SeriesCoefficient[SeriesCoefficient[-3*s*u*(s + u), {u,0,2}], {s,0,1}]//FullSimplify]

Print["su^2 coeff = ", SeriesCoefficient[SeriesCoefficient[-3*s*u*(s + u), {u,0,0}], {s,0,3}]//FullSimplify]

