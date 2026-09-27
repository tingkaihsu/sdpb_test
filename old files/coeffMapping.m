(* ::Package:: *)

fwdM[s_,u_,m_]:= g0 + g2*((s-2m^2)^2+(u-2m^2)^2)+d2*t^2+g3*((s-2m^2)^3+(u-2m^2)^3)+d3*t^3/.{t->4m^2-s-u};

Print[fwdM[s,u,m]//Simplify]

stuM[s_,u_,m_]:=G0+G2*((s-4m^2/3)^2+(u-4m^2/3)^2+(t-4m^2/3)^2)+G3((t-4m^2/3)^3+(s-4m^2/3)^3+(u-4m^2/3)^3)/.{t->4m^2-s-u};

Print[stuM[s,u,m]//Simplify]

Print["G2 = ", SeriesCoefficient[stuM[s,u,m],{s,4m^2/3,2},{u,4m^2/3,0}]]
Print["G2 = ", SeriesCoefficient[stuM[s,u,m],{s,4m^2/3,1},{u,4m^2/3,1}]]
Print["G2 = ", SeriesCoefficient[stuM[s,u,m],{s,4m^2/3,0},{u,4m^2/3,2}]]

Print["g[2,0] = ", SeriesCoefficient[stuM[s,u,m],{s,2m^2,2},{u,2m^2,0}]]
Print["g[1,1] = ", SeriesCoefficient[stuM[s,u,m],{s,2m^2,1},{u,2m^2,1}]]
Print["g[0,2] = ", SeriesCoefficient[stuM[s,u,m],{s,2m^2,0},{u,2m^2,2}]]

Print["g[1,1]-2g[2,0] = ", SeriesCoefficient[stuM[s,u,m],{s,2m^2,1},{u,2m^2,1}]-2*SeriesCoefficient[stuM[s,u,m],{s,2m^2,2},{u,2m^2,0}]//Simplify]



