(* ::Package:: *)

gABABCoeff[a_, b_] := gABAB @@ Sort[{a, b}];

validTriples[Nmax_Integer] :=
    Flatten[Table[{a, b}, {a, 0, Nmax}, {b, 0, Nmax - a}], 1];
    
MABAB[s_, t_, mA_, Nmax_Integer] := Module[{u},
	u = 2mA^2-s-t;
    -gABB^2(1/s+1/u)-gAAA*gABB*(1/(t-mA^2))+
    Total[
        Function[{ab},
            gABABCoeff[ab[[1]], ab[[2]]]
            * (s-mA^2)^ab[[1]]
            * (u-mA^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] ];
    

gAABBCoeff[a_,b_] := gABAB @@ Sort[{a,b}];

MAABB[s_,t_,m_,Nmax_Integer] :=Module[{u},
	u = 2m^2-s-t;
    -gAAA*gABB(1/(s-m^2))-gABB^2(1/t+1/u)+
    Total[
        Function[{ab},
            gABABCoeff[ab[[1]], ab[[2]]]
            * (t-m^2)^ab[[1]]
            * (u-m^2)^ab[[2]]
        ] /@ validTriples[Nmax]
    ] ];

(* Series[-Residue[MABAB[s,t,m,10]/((s-m^2)*((s-m^2)*(s-m^2+t))^(k/2))/.{k->0},{s,Infinity}],{t,0,0}] *)

MABAB[s,t,m,10]-MAABB[t,s,m,10]//FullSimplify

(* test Kernel *)

KerLow[s0_,m_,k_,q_,Nmax_Integer]:= Residue[Residue[MABAB[s,t,m,Nmax]/((s-s0)^(k-q+1)t^(q+1)),{s,Infinity}],{t,0}]-Residue[Residue[MAABB[s,t,m,Nmax]/((t-s0)^(k-q+1)s^(q+1)),{s,Infinity}],{t,s0}];


(* sum rule *)

(* partial waves *)
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
(* Matched AA/BB-channel convention: Sqrt[PhiAA PhiBB], with PhiAA = (x-4m^2)^(7/2)/Sqrt[x] (AAAA file)
   and PhiBB = x^3 (BBBB file). Keep these three consistent; see the note above. *)
phiAABB[s_,m_]:= s^(5/4) (s - 4 m^2)^(7/4);

(* without phase-volume factor *)
KerABAB[x_,J_,k_,q_,m_]:=Residue[(1/((x+m^2)^(k-q+1)*t^(q+1))-1/((2m^2-x-t+m^2)^(k-q+1)*t^(q+1)))*PartialWaveD[10,J,1 + 2 x t/(x - m^2)^2],{t,0}];


(* term 2: s' on the neutral cut ... *)
KerAABB1[x_,J_,k_,q_,m_]:=Residue[(1/((x)^(q+1)*(t+m^2)^(k-q+1)))*PartialWaveD[10,J,(2 t + x - 2 m^2)/Sqrt[x (x - 4 m^2)]],{t,-m^2}];

(* term 3: the second term on the u-channel cut *)
KerAABB2[x_,J_,k_,q_,m_]:=Residue[-(1/((2m^2-x-t)^(q+1)*(t+m^2)^(k-q+1)))*(-1)^J*PartialWaveD[10,J,-1-(2 (2 m^2-t-x) x)/(-m^2+x)^2],{t,-m^2}];


NABAB[x_,m_,J_,k_, q_] := KerABAB[x, J, k, q, m] - KerAABB2[x, J, k, q, m];
NAABB[x_,m_,J_,k_, q_] := -KerAABB1[x, J, k, q, m];


(* ---------- index list: k,q ---------- *)
kMax = 8;
nullIndex = Flatten[Table[{k, q}, {k, 0, kMax}, {q, 0, k}], 1];
(* literal pseudo code would be {q, 0, Floor[k/2]}: adds 10 identically-zero rows *)
Length[nullIndex]

(* ---------- the double loop ---------- *)
LaunchKernels[];

nullABABx = ParallelMap[Together@NABAB[z, mA, J, #[[1]], #[[2]]] &, nullIndex];
nullAABBx = ParallelMap[Together@NAABB[z, mA, J, #[[1]], #[[2]]] &, nullIndex];

ABABlowPole[s_, t_, mA_] := -gABB^2(1/s+ 1/(2mA^2-s-t) )-gAAA*gABB*(1/(t-mA^2));

AABBlowPole[s_, t_, mA_] := -gAAA*gABB(1/(s-m^2))-gABB^2(1/t+1/(2mA^2-s-t) );

DblCtrABAB[k_Integer, q_Integer] := Residue[Residue[
    1/( (s+mA^2) * t ) * ABABlowPole[s,t,mA]/( (s+mA^2)^(k - q) * (t)^q ),
    {s, Infinity}
], {t, 0} ];

(* with additional minus sign *)
DblCtrAABB[k_Integer, q_Integer] := Residue[Residue[
    -1/( (t+mA^2) * s ) * AABBlowPole[s,t,mA]/( (t+mA^2)^(k - q) * (s)^q ),
    {s, Infinity}
], {t, -mA^2}]

nullABABa = ParallelMap[Together@DblCtrABAB[#[[1]], #[[2]]] &, nullIndex];

nullAABBa = ParallelMap[Together@DblCtrAABB[#[[1]], #[[2]]] &, nullIndex];


Length[nullAABBx]


ToString[nullAABBx,InputForm,PageWidth->Infinity]


ToString[nullAABBa/.{gABB->1},InputForm,PageWidth->Infinity]


ToString[nullABABa/.{gABB->1,gAAA->1},InputForm,PageWidth->Infinity]


(* ======================================================================== *)
(* EXPORT.  Gegenbauer-form kernels for the PMP generators (test20.m)       *)
(*                                                                          *)
(* Writes kernels_AB.m with                                                 *)
(*   nGeg[n, l, z]          normalized Gegenbauer, P_n^(l)(1) = 1, 0 if n<0 *)
(*   NABAB[n_,z_,J_,m_]     = nullABABx[[n+1]]   (same format as before)    *)
(*   NAABB[n_,z_,J_,mA_]    = nullAABBx[[n+1]]                              *)
(* Every Hypergeometric2F1[k-J, J+7+k, 4+k, w] equals nGeg[J-k, 7/2+k, 1-2w]*)
(* (P_J^(k) = P_J^(k)(1) * normalized C_{J-k}^{7/2+k}).  Each distinct 2F1  *)
(* is computed ONCE per call through With[] (NABAB has 9814 calls but only  *)
(* 9 distinct values), and GegenbauerC keeps full precision at J ~ 5*10^4,  *)
(* where the numerical Hypergeometric2F1 loses all digits near threshold.   *)
(* ======================================================================== *)

nGegDefinition =
  "nGeg[n_Integer, l_, z_] := If[n < 0, 0, GegenbauerC[n, l, z] n!/Pochhammer[2 l, n]];";

kernelDefinitionString[lhs_String, list_List] := Module[{hyps, bad, names, body, locals},
  hyps = Union @ Cases[list, _Hypergeometric2F1, {0, Infinity}];
  (* only the family 2F1[k-J, J+7+k, 4+k, w] is allowed *)
  bad = Select[hyps,
    Expand[#[[1]] + #[[2]] - 2 #[[3]] + 1] =!= 0 || !IntegerQ[Expand[#[[1]] + J]] &];
  If[bad =!= {}, Print["kernelDefinitionString: unexpected 2F1 ", bad]; Abort[]];
  names  = Table["gg" <> ToString[i], {i, Length[hyps]}];
  ClearAll @@ names;          (* placeholders must be free symbols *)
  body   = list /. Thread[hyps -> (Symbol /@ names)];
  locals = MapThread[
    #1 <> " = nGeg[" <> ToString[-#2[[1]], InputForm] <> ", " <>
      ToString[#2[[3]] - 1/2, InputForm] <> ", " <>
      ToString[1 - 2 #2[[4]], InputForm, PageWidth -> Infinity] <> "]" &,
    {names, hyps}];
  lhs <> " := With[{" <> StringRiffle[locals, ", "] <> "}, " <>
    ToString[body, InputForm, PageWidth -> Infinity] <> "][[n+1]];"
];

exportKernels[file_String] := Export[file,
  StringRiffle[{
    "(* kernels_AB.m -- generated by ABAB_check.m, kMax = " <> ToString[kMax] <>
      ", " <> DateString[] <> ".  Do not edit by hand. *)",
    nGegDefinition,
    kernelDefinitionString["NABAB[n_,z_,J_,mA_]",  nullABABx],
    kernelDefinitionString["NAABB[n_,z_,J_,mA_]", nullAABBx]
  }, "\n\n"], "Text"];

kernelsFile = FileNameJoin[{
  Which[
    StringQ[$InputFileName] && StringLength[$InputFileName] > 0,     (* wolframscript / Get *)
      DirectoryName[ExpandFileName[$InputFileName]],
    $FrontEnd =!= Null,                                                (* run from the front end *)
      Quiet @ Check[NotebookDirectory[], Directory[]],
    True, Directory[]],
  "kernels_AB.m"}];
exportKernels[kernelsFile];
Print["wrote ", kernelsFile];
