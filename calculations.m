result = 1/3* ( Integrate[(-2*t+2)*Exp[-I*m*2Pi/3*t], {t,0,1}] + Integrate[(t-1)*Exp[-I*m*2Pi/3*t], {t,1,3}] )

Print[ FullSimplify[result] ]

(3 - 9*E^((4*I)/3*m*Pi) + (4*I)*m*Pi + E^((2*I)*m*Pi)*(6 - (4*I)*m*Pi))/(4*E^((2*I)*m*Pi)*m^2*Pi^2)
