(* ::Package:: *)

gABABCoeff[a_, b_] := gABAB @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
      Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];
    
MABAB[s_, t_, mA_, Nmax_Integer] := Module[{u},
	  u = 2mA^2-s-t;
      Total[
            Function[{ab},
                  gABABCoeff[ab[[1]], ab[[2]]]
                   * (s - mA^2)^ab[[1]]
                   * (u - mA^2)^ab[[2]]
              ] /@ validTriples[Nmax]
        ]];

KerCoeff[s_,t_,m_,k_, Nmax_Integer,s1_,s2_] := MABAB[s, t, m, Nmax]/((s-s1)((s-s1)*(s-s2))^(k/2))


s1 = m^2;
s2 = m^2-t;

(* Note the additional minus sign *)
-SeriesCoefficient[Residue[KerCoeff[s,t,m,0,10,s1,s2],{s,Infinity}],{t,0,0}]


PartialWaveD[d_, J_, z_] := Hypergeometric2F1[-J, J + d - 3, (d - 2)/2, (1 - z)/2];
(* AB-channel weight convention (positivity-only bounds).
   The code multiplies by the AB phase space PhiAB = (x - m^2)^7/x^4 (d = 10, m_B = 0),
   so the AB density it implicitly uses is rho~_J = n_J rho_J / PhiAB^2 >= 0,
   where Im M_AB = Sum_J n_J rho_J/PhiAB P_J is the physical normalization.
   This is equivalent to the physical setup ONLY if every AB-channel row uses wAB:
   SumABAB (gAB00.m), KerABAB and KerAABB2 (ABAB_check.m), and the NABAB rows in the SDPB files.
   Switch to 1/PhiAB everywhere before imposing unitarity upper bounds (rho_J <= 2, |S_J| <= 1)
   or comparing spectral densities to explicit UV models. *)
wAB[x_, m_] := (x - m^2)^7/x^4;

SumABAB[x_,t_,J_,k_,m_,x1_,x2_]:= (1/(x-x1)*1/((x-x1)*(x-x2))^(k/2)-1/((2m^2-x-t-x1)*((2m^2-x-t-x1)*(2m^2-x-t-x2))^(k/2)))*wAB[x, m]*PartialWaveD[10,J,1 + 2 x t/(x - m^2)^2];




SeriesCoefficient[SumABAB[x,t,J,0,m,s1,s2],{t,0,0}]//FullSimplify
