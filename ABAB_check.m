(* ::Package:: *)

ParentAmp[sh_, th_] := (-Gamma[-th]*Gamma[-uh])/Gamma[1+sh]-(Gamma[-sh]*Gamma[-uh])/Gamma[1+th]-(Gamma[-sh]*Gamma[-th])/Gamma[1+uh]/.{uh -> -sh-th};
  
TestABAB[s_, t_, m_] :=
  ParentAmp[
    s - m^2,
    t
  ];

TestAABB[s_, t_, m_] :=
  ParentAmp[
    s,
    t - m^2
  ];

(* Test s-t crossing *)
TestABAB[s,t,m]-TestAABB[t,s,m]


KerABAB[s_,t_,m_,k_,q_]:=1/(s^(k-q+1)t^(q+1))*TestABAB[s+m^2,t,m];

KerAABB[s_,t_,m_,k_,q_]:=1/(s^(q+1)t^(k-q+1))*TestAABB[s,m^2+t,m];


DblCtrCheck[m_,k_,q_]:=Residue[Residue[KerABAB[s,t,m,k,q]-KerAABB[s,t,m,k,q],{s,0}],{t,0}]


DblCtrCheck[m,1,0]
