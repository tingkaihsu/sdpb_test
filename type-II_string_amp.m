(* ============================================================
   Type-II string amplitude: KK reduction 5d -> 4d
   5d dilaton scattering -> 4d massive scalar amplitude
   ============================================================ *)

(* 5d type-II closed string amplitude (Virasoro-Shapiro-like).
   u is a free symbol here; the Mandelstam constraint u=-s-t
   is enforced by the replacement rule. *)
A4[s_,t_] := (s^2+t^2+u^2)^2 *
              (Gamma[-s/4] Gamma[-t/4] Gamma[-u/4]) /
              (Gamma[1+s/4] Gamma[1+t/4] Gamma[1+u/4])/. {u -> -s-t};

Print[ A4[s,t] // FullSimplify ];

(* s-channel *)
s5d[s_, m_] := s - (m+m)^2;
t5d[t_, m_] := t;
u5d[u_, m_] := u;

(* Verify: shifted sum is 4m^2 on the Mandelstam constraint surface *)
Print[ s5d[s,m] + t5d[t,m] + u5d[u,m] /. {s+t+u -> 4m^2} ];

Print["A4[s5d[s,m], t] = ", A4[s5d[s,m], t]//FullSimplify ]

M4[s_, t_, m_] := ( (s-4m^2)^2+t^2+u^2 )^2 *
                  ( Gamma[-(s-4m^2)/4] Gamma[-t/4] Gamma[-u/4] ) /
                  ( Gamma[1+(s-4m^2)/4] Gamma[1+t/4] Gamma[1+u/4] ) 
                  + ( s^2+(t-4m^2)^2+u^2 )^2 *
                  ( Gamma[-s/4] Gamma[-(t-4m^2)/4] Gamma[-u/4] ) /
                  ( Gamma[1+s/4] Gamma[1+(t-4m^2)/4] Gamma[1+u/4] )
                  + ( s^2+t^2+(u-4m^2)^2 )^2 *
                  ( Gamma[-s/4] Gamma[-t/4] Gamma[-(u-4m^2)/4] ) /
                  ( Gamma[1+s/4] Gamma[1+t/4] Gamma[1+(u-4m^2)/4] );
(* ============================================================
   RESIDUE COMPUTATION \[LongDash] CORRECTED

   The original code called:
     Residue[M4[s,t]/((s-2m^2)^2 (s-(2m^2-t))), {s, 2m^2}]
     Residue[M4[s,t]/((s-2m^2)^2 (s-(2m^2-t))), {s, 2m^2-t}]

   WHY THIS IS EXTREMELY SLOW:
     Residue[] internally calls Series[] to Laurent-expand the
     full expression around the given point. Here the expression
     contains Gamma[-s/4], Gamma[-t/4], and Gamma[-u(s)/4].
     Expanding these symbolically around s=2m^2 (with generic m)
     produces PolyGamma series \[LongDash] nine such expansions across three
     terms. Mathematica cannot simplify this in finite time for
     generic m.

   THE KEY OBSERVATION:
     M4[s,t] is ANALYTIC at s=2m^2 for generic m. The Gamma
     function Gamma[-s/4] has poles only at s=0,4,8,...  (i.e.
     when -s/4 is a non-positive integer). For generic m, 2m^2
     is none of these. Therefore the poles of the integrand
     g(s) = M4[s,t] / ((s-2m^2)^2 * (s-(2m^2-t)))
     come ENTIRELY from the explicit polynomial denominator.

   ANALYTIC RESIDUE FORMULAS:
     Let  a = 2m^2,   b = 2m^2-t,   so  a-b = t.

     Double pole at s=a:
       Res[g, a] = d/ds [M4[s,t]/(s-b)] |_{s=a}
                 = (M4'(a,t)*(a-b) - M4(a,t)) / (a-b)^2
                 = (M4'(2m^2,t)*t - M4[2m^2,t]) / t^2

     Simple pole at s=b:
       Res[g, b] = M4[b,t] / (b-a)^2
                 = M4[2m^2-t, t] / t^2

   This requires only three substitution/differentiation steps \[LongDash]
   no Series expansion of Gamma functions whatsoever.
   ============================================================ *)

(* Step 1: Evaluate M4 at the three required points *)
M4at    = M4[2m^2, t, m]       // FullSimplify;
M4deriv = D[M4[s,t, m], s] /. s -> 2m^2 // FullSimplify;
M4shft  = M4[2m^2-t, t, m]     // FullSimplify;

(* Step 2: Assemble the residue from the analytic formula *)
residue = (M4deriv * t - M4at + M4shft) / t^2 // FullSimplify;

(* Step 3: Series coefficients in t *)
Print[SeriesCoefficient[residue, {t, 0, 0}] ]
Print[SeriesCoefficient[residue, {t, 0, 1}] ]