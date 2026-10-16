# AA → AA null constraints: derivation and explicit block matrices

Conventions match `test20.m` (option R, residue units, test16 treatment of light states).
Every kernel below was checked symbolically against `DblCtrSum[x,k,q]` of `AAAA_sum_rules.m`
divided by $\phi_{AA}$ (exact agreement at $J=0,2,4,6$ and several $x$).

---

## 1. Setup

Kinematics and shorthand:

$$
s+t+u = 4m^2,\qquad m^2\equiv m_A^2,\qquad y\equiv x-4m^2,\qquad \mathcal J\equiv J(J+7).
$$

A neutral state of mass$^2$ $x$ and even spin $J$ contributes, near its $s$-channel pole,

$$
M(s,t)\;\simeq\;-\,R\,\frac{P_J(z)}{s-x},\qquad
z = 1+\frac{2t}{x-4m^2},\qquad
P_J(z)={}_2F_1\!\Big(-J,\,J+7;\,4;\,\tfrac{1-z}{2}\Big),\quad P_J(1)=1,
$$

with residue $R\ge 0$. The spectrum is

| state | $x$ | $J$ | residue $R$ |
|---|---|---|---|
| $A$ (light) | $m^2$ | $0$ | $g_{AAA}^2$ |
| $H$ | $m_2^2=1$ | $2$ | $\phi_{AA}(x)\,\tilde\rho_{AA}$ |
| continuum | $x\ge m_{\rm gap}^2$ | even | $\phi_{AA}(x)\,\tilde\rho_{AA}$ |

$$
\phi_{AA}(x)=\frac{(x-4m^2)^{7/2}}{\sqrt{x}} .
$$

The low-energy amplitude is a pure crossing-symmetric polynomial (no IR poles); $A$ is an ordinary state.

**Assumptions.** (A1) crossing $M(s,t)=M(t,s)=M(u,t)$; (A2) narrow resonances;
(A3) $M$ analytic near $(s,t)=(0,0)$, i.e. $4m^2<\min(m_2^2,m_{\rm gap}^2)$ and $4m^2\neq m_2^2$;
(A4) Regge: $|M|\lesssim |s|^{\alpha_N(t)}$ with $\alpha_N(0)<0$
(the trajectory through $A$ and $H$ gives $\alpha_N(0)=-2m^2/(m_2^2-m^2)=-0.5$).

---

## 2. Fixed-$t$ dispersion relation

Expand at $(s,t)=(0,0)$, $u_0=4m^2$: $M(s,t)=\sum_n a_n(t)\,s^n$.
Deforming $\oint ds\,M/s^{n+1}$ to infinity (arc $\to 0$ by A4 for all $n\ge 0$):

- $s$-channel pole: $\;-R P/(s-x)=R P/(x-s)\;\Rightarrow\; R P/x^{n+1}$
- $u$-channel pole: $\;u=x \Leftrightarrow s=4m^2-x-t\;\Rightarrow\; -R P/(4m^2-x-t)^{n+1}$

$$
\boxed{\,a_n(t)=\sum_{\rm states} R\;F_n(x,J;t),\qquad
F_n = P_J(z)\left[\frac{1}{x^{n+1}}-\frac{1}{(4m^2-x-t)^{n+1}}\right]\,}
$$

(identical particles, even $J$: the $u$-channel angle gives the same $P_J(z)$).

---

## 3. Null constraints from $s\leftrightarrow t$

With $M=\sum c_{n,q}\,s^n t^q$, $c_{n,q}=[t^q]\,a_n(t)$, crossing $c_{n,q}=c_{q,n}$ gives, for $k=n+q$:

$$
\boxed{\,\sum_{\rm states} R\;K_{k,q}(x,J)=0,\qquad
K_{k,q}=[t^q]\,F_{k-q}-[t^{k-q}]\,F_q\,}
$$

$K_{k,k-q}=-K_{k,q}$ and $K_{k,k/2}=0$, so only $0\le q\le\lfloor (k-1)/2\rfloor$ is kept.
This is exactly `DblCtrSum[x,k,q]`$/\phi_{AA}$.

**Closed form.** With $n=k-q$, $u_0=4m^2-x=-y$,

$$
K_{k,q}=\frac{p_q}{x^{n+1}}-\frac{p_n}{x^{q+1}}
-\sum_{j=0}^{q}\binom{n+j}{j}\frac{p_{q-j}}{u_0^{\,n+1+j}}
+\sum_{j=0}^{n}\binom{q+j}{j}\frac{p_{n-j}}{u_0^{\,q+1+j}},
$$

where $p_j=[t^j]\,P_J(z)$:

$$
p_j=\frac{(-J)_j\,(J+7)_j}{(4)_j\,j!}\left(-\frac{1}{y}\right)^{j},\qquad
(-J)_j(J+7)_j=(-1)^j\prod_{i=0}^{j-1}\big(\mathcal J-i(i+7)\big),
$$

$$
p_0=1,\quad p_1=\frac{\mathcal J}{4y},\quad p_2=\frac{\mathcal J(\mathcal J-8)}{40\,y^2},\quad
p_3=\frac{\mathcal J(\mathcal J-8)(\mathcal J-18)}{720\,y^3},\;\dots
$$

---

## 4. The first six kernels ($y=x-4m^2$, $\mathcal J=J(J+7)$, $x+y=2(x-2m^2)$, $x-y=4m^2$)

$$
K_{1,0}=\frac{1}{x^2}-\frac{\mathcal J\,(x+y)}{4\,x\,y^2}
$$

$$
K_{2,0}=\frac{1}{x^3}-\frac{\mathcal J\big[\mathcal J(x+y)-18x-8y\big]}{40\,x\,y^3}
$$

$$
K_{3,0}=\frac{1}{x^4}-\frac{\mathcal J\big[\mathcal J^2(x+y)-\mathcal J(44x+26y)+468x+144y\big]}{720\,x\,y^4}
$$

$$
K_{3,1}=\frac{4m^2\,\mathcal J\big[(\mathcal J-18)\,x(x+y)-10\,y^2\big]}{40\,x^3\,y^4}
$$

$$
K_{4,0}=\frac{1}{x^5}-\frac{\mathcal J\big[\mathcal J^3(x+y)-\mathcal J^2(84x+56y)+\mathcal J(2156x+924y)-17424x-4320y\big]}{20160\,x\,y^5}
$$

$$
K_{4,1}=\frac{\mathcal J\big[\mathcal J^2x^2(x^2-y^2)-\mathcal J(62x^4-26x^2y^2)+792x^4-144x^2y^2+180y^4\big]}{720\,x^4\,y^5}
$$

**Large $J$ (fixed $x$).** $K_{k,0}\simeq-\dfrac{\mathcal J^{k}\,(x+y)}{(4)_k\,k!\;x\,y^{k+1}}$, i.e. growth $J^{2k}$ with a negative sign;
$K_{k,q}$ grows like $J^{2(k-q)}$ at most. This is what the large-$J$ / impact blocks of test16 must control.

**Forward weight (the $g_{0AA}$ row, $n=q=0$):**

$$
w_{AA}(x)=\frac1x+\frac{1}{x-4m^2}=\frac{x+y}{x\,y}>0\quad(x>4m^2).
$$

---

## 5. The light state $A$ ($x=m^2$, $J=0$, $y=-3m^2$)

All binomial terms cancel between the two halves of $K$:

$$
K_{k,q}(m^2,0)=\frac{\delta_{q0}-\delta_{qk}}{m^{2(k+1)}}\;\xrightarrow{\;q\le (k-1)/2\;}\;\frac{\delta_{q0}}{m^{2(k+1)}},
\qquad w_{AA}(m^2)=\frac{2}{3m^2}.
$$

| row | $g_{0AA}$ | $N_{1,0}$ | $N_{2,0}$ | $N_{3,0}$ | $N_{3,1}$ | $N_{4,0}$ | $N_{4,1}$ |
|---|---|---|---|---|---|---|---|
| general | $\frac{2}{3m^2}$ | $\frac1{m^4}$ | $\frac1{m^6}$ | $\frac1{m^8}$ | $0$ | $\frac1{m^{10}}$ | $0$ |
| $m^2=\tfrac15$ | $\tfrac{10}{3}$ | $25$ | $125$ | $625$ | $0$ | $3125$ | $0$ |

These equal the exact double contour of the physical-sign pole $-g^2\big[\tfrac1{s-m^2}+\tfrac1{t-m^2}+\tfrac1{u-m^2}\big]$, moved to the states side.

---

## 6. Explicit block matrices (coordinates restricted to the AA part)

Coordinate basis (the BBBB and cross rows are elided as "$\cdots$"):

$$
\vec e=\big(\,g_{0AA},\;g_{0BB},\;X_H\;\big|\;N_{1,0},\,N_{2,0},\,N_{3,0},\,N_{3,1},\,N_{4,0},\,N_{4,1},\dots\big|\;\cdots\big).
$$

### 6.1 Heavy neutral state (`PolyH`, `PolyN`), basis $\{AA,BB\}$

Each element carries one phase-space factor on **all** its components
($\phi_X=\sqrt{\phi_{AA}\phi_{BB}}$, $\xi=$ `xhSign` for $H$, $0$ for the continuum):

$$
\mathcal M^{(x,J)}=
\begin{pmatrix}
\phi_{AA}(x)\,\vec v_{AA} & \phi_X(x)\,\vec v_X\\[2pt]
\phi_X(x)\,\vec v_X & \phi_{BB}(x)\,\vec v_{BB}
\end{pmatrix},
$$

$$
\begin{aligned}
\vec v_{AA}&=\big(\,w_{AA},\;0,\;0\;\big|\;K_{1,0},\,K_{2,0},\,K_{3,0},\,K_{3,1},\,K_{4,0},\,K_{4,1},\dots\big|\;0\;\big),\\
\vec v_{BB}&=\big(\,0,\;w_{BB},\;0\;\big|\;0,\,0,\,0,\,0,\,0,\,0,\dots\big|\;\text{BBBB rows}\;\big),\\
\vec v_{X}&=\big(\,0,\;0,\;\tfrac{\xi}{2}\;\big|\;0,\,0,\,0,\,0,\,0,\,0,\dots\big|\;-\tfrac12\,\text{cross rows}\;\big).
\end{aligned}
$$

Row by row, an AA→AA null constraint $N_{k,q}$ sees only the $AA$ element:

$$
\mathcal M^{(x,J)}_{N_{k,q}}=
\begin{pmatrix}\phi_{AA}(x)\,K_{k,q}(x,J) & 0\\ 0 & 0\end{pmatrix},
\qquad
\mathcal M^{(x,J)}_{g_{0AA}}=
\begin{pmatrix}\phi_{AA}(x)\,w_{AA}(x) & 0\\ 0 & 0\end{pmatrix}.
$$

Explicitly, e.g. for the first row,

$$
\phi_{AA}K_{1,0}=\frac{(x-4m^2)^{7/2}}{\sqrt x}\left[\frac1{x^2}-\frac{J(J+7)\,(x-2m^2)}{2\,x\,(x-4m^2)^2}\right],
\qquad
\phi_{AA}w_{AA}=\frac{2(x-2m^2)(x-4m^2)^{5/2}}{x^{3/2}} .
$$

At $H$ ($x=1$, $J=2$, $\mathcal J=18$, $m^2=\tfrac15$): $\phi_{AA}=0.2^{7/2}\approx3.58\times10^{-3}$,
$\phi_X=0.2^{7/4}\approx5.98\times10^{-2}$, $\phi_{BB}=1$.

### 6.2 Light block (`PolyA`), basis $(g_{AAA},g_{ABB})$, weight 1

$$
\mathcal M^{(A)}=
\begin{pmatrix}
\vec a & \vec b\\ \vec b & \vec c
\end{pmatrix},
\qquad
\vec a=\Big(\tfrac{2}{3m^2},\,0,\,0\;\Big|\;\tfrac1{m^4},\,\tfrac1{m^6},\,\tfrac1{m^8},\,0,\,\tfrac1{m^{10}},\,0,\dots\Big|\;0\Big),
$$

while $\vec b$ (cross rows) and $\vec c$ ($g_{0BB}$, BBBB and cross rows) are zero on all AA→AA rows.
Row by row:

$$
\mathcal M^{(A)}_{N_{k,q}}=\begin{pmatrix}\delta_{q0}\,m^{-2(k+1)} & 0\\ 0&0\end{pmatrix},\qquad
\mathcal M^{(A)}_{g_{0AA}}=\begin{pmatrix}\tfrac{2}{3m^2} & 0\\ 0&0\end{pmatrix}.
$$

No $\phi$ here: $\phi_{AA}(m^2)$ is complex ($A$ lies below the $AA$ threshold); only the residue exists.

### 6.3 The resulting equations

$$
N_{k,q}:\quad
g_{AAA}^2\,\frac{\delta_{q0}}{m^{2(k+1)}}
+\phi_{AA}(m_2^2)\,\tilde\rho^{H}_{AA}\,K_{k,q}(m_2^2,2)
+\sum_{x\ge m_{\rm gap}^2}\sum_{J\ \rm even}\phi_{AA}(x)\,\tilde\rho_{AA}(x,J)\,K_{k,q}(x,J)=0,
$$

$$
g_{0AA}=\frac{2\,g_{AAA}^2}{3m^2}+\sum_{\rm heavy}\phi_{AA}(x)\,\tilde\rho_{AA}(x,J)\,w_{AA}(x).
$$

Every term is (residue) $\times$ (stripped kernel), so the heavy $\phi$-units and the light residue units add consistently.

---

## 7. Checks performed

1. $\phi_{AA}K_{k,q}$ = `DblCtrSum[x,k,q]` (AAAA_sum_rules.m) exactly, for all six labels at $(x,J)=(3,2),(7/2,4),(5,6),(1,0)$, $m^2=1/5$.
2. Light-state values reproduce $\delta_{q0}/m^{2(k+1)}$ and $w_{AA}(m^2)=2/(3m^2)$.
3. All 42 labels with $q\le\lfloor(k-1)/2\rfloor$, $k\le 12$, are linearly independent (exact rank, $m^2=1/5$).
