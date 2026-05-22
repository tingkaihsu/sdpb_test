(* ::Package:: *)

(* ::Package:: *)
(**)


(* ============================================================
   Type-II string amplitude: KK reduction 5d -> 4d
   5d dilaton scattering -> 4d massive scalar amplitude
   ============================================================ *)

(* 5d type-II closed string amplitude (Virasoro-Shapiro-like).
   u is a free symbol here; the Mandelstam constraint u=-s-t
   is enforced by the replacement rule. *)
A4[s_,t_] := (s^2+t^2+u^2)^2 *
              (Gamma[-s/4] Gamma[-t/4] Gamma[-u/4]) /
              (Gamma[1+s/4] Gamma[1+t/4] Gamma[1+u/4]) /. {u -> -s-t};

A4[s,t] // FullSimplify


mstr = n/R;

(* Mandelstam variables shifted by KK mass in each channel.
   The KK momentum along the compact direction contributes
   (m+m)^2 = 4m^2 to the s-channel invariant. *)
hs[s_, m_] := s + (m+m)^2;   (* = s + 4m^2 *)
ht[t_, m_] := t;
hu[u_, m_] := u;

(* Verify: shifted sum is 4m^2 on the Mandelstam constraint surface *)
hs[s,m] + ht[t,m] + hu[u,m] /. {s+t+u -> 0}


(* s-channel shifted amplitude (as a sanity check) *)
A4[hs[s,m], t]


(* ============================================================
   stu-symmetric 4d amplitude.

   Three terms, one per channel:
     term 1: s-channel  -- shift s -> s+4m^2, u_eff = -(s+4m^2)-t  (inside A4)
     term 2: t-channel  -- shift t -> t+4m^2, u_eff = -s-(t+4m^2)  (inside A4)
     term 3: u-channel  -- u -> 4m^2-s-t     (written out explicitly)

   BUG FIX (minor, operator-precedence / clarity):
   The original code wrote:
     A4[...] + A4[...] + (Gamma-expr) /. {u -> 4m^2-s-t}
   Mathematica's /. has LOWER precedence than +, so this is
   equivalent to (A4[...]+A4[...]+Gamma-expr) /. {u->4m^2-s-t}.
   In practice this is harmless (u is already gone from the first
   two terms after A4 evaluates), but it is fragile and misleading.
   Fix: wrap the third term in explicit parentheses.
   ============================================================ *)
M4[s_,t_] :=
  A4[s+4m^2, t] +
  A4[s, t+4m^2] +
  ((s^2+t^2+u^2)^2 *
   (Gamma[-s/4] Gamma[-t/4] Gamma[-u/4]) /
   (Gamma[1+s/4] Gamma[1+t/4] Gamma[1+u/4]) /. {u -> 4m^2-s-t});


(* Symmetry checks *)
Print["s-t symmetry test of A4[s,t] = ",
      A4[s+4m^2,t] - A4[t+4m^2,s] // FullSimplify]
Print["s-t symmetry test of M4[s,t] = ",
      M4[s,t] - M4[t,s]]


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
M4at[m_, t_]    := M4[2m^2, t];
M4deriv[m_, t_] := D[M4[s,t], s] /. s -> 2m^2;
M4shft[m_,t_]   := M4[2m^2-t, t];

(* Step 2: Assemble the residue from the analytic formula *)
residue[m_, t_] := (M4deriv[m,t] * t - M4at[m,t] + M4shft[m,t]) / t^2;

(* Step 3: Series coefficients in t *)
g2[m_] := SeriesCoefficient[residue[m,t], {t, 0, 0}];
g3[m_] := SeriesCoefficient[residue[m,t], {t, 0, 1}];



Plot[g3[m]/g2[m],{m,0,4/10}]
