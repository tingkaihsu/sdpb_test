"""
Sampling grids of test20.m (b_max = 100 everywhere) and test20v3.m
(large-J tail extended to b = 200), mAsq = 1/5, drawn in the (M^2, J) plane.

Formulas are transcribed from test20v3.m (BLOCK 1, suppMasses, PMP2SDP):
  base     M^2 = gap / (1 - sin^2(phi/2)),  phi_k = pi/2 + pi/2 cos((k+1/2) pi/150)
  supp     J >= 40, theta_j = j dPsi/(J+7/2), j = 1..K,
           K = floor(min((J+7/2) theta_s(gap), m_A bmax(J)) / dPsi)
  impact   b_l = l dPsi/m_A, l = 0..floor(m_A b_imp/dPsi), plus b_imp
test20: bmax(J) = 100, b_imp = 100.   test20v3: bmax(J) = 200 if J >= 400 else 100, b_imp = 200.

Run:  python plot_sampling.py      (writes fig*.png next to this file)
"""
import os
from math import pi, sin, cos, asin, acos, sqrt, floor
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
mAsq = 1 / 5
mA = sqrt(mAsq)
dPsi = pi / 4
npts = 150
GAP = {"N": 166 / 100, "C": 1.0}
SPINS = {
    "N": list(range(0, 201, 2)) + [400, 1000, 10000, 50000],
    "C": list(range(0, 201)) + [x for y in (400, 1000, 10000, 50000) for x in (y, y + 1)],
}
NAME = {"N": "neutral (J even, $M^2\\geq 1.66$)", "C": "charged (all J, $M^2\\geq 1$)"}


def theta(s, M2):
    return asin(2 * sqrt(mAsq / M2)) if s == "N" else acos(1 - 2 * mAsq * (M2 + mAsq) / (M2 - mAsq) ** 2)


def mass_of_theta(s, th):
    if s == "N":
        return 4 * mAsq / sin(th) ** 2
    c = 1 - cos(th)
    return mAsq * ((c + 1) + sqrt(4 * c + 1)) / c


def base(s):
    phi = pi / 2 + pi / 2 * np.cos((np.arange(npts) + 0.5) * pi / npts)
    return GAP[s] / (1 - np.sin(phi / 2) ** 2)


def supp(s, J, bmax):
    if J < 40:
        return []
    nn = J + 3.5
    K = floor(min(nn * theta(s, GAP[s]), mA * bmax) / dPsi)
    return [mass_of_theta(s, j * dPsi / nn) for j in range(1, K + 1)]


def bmax_v20(J):
    return 100


def bmax_v3(J):
    return 200 if J >= 400 else 100


def impact_nodes(bimp):
    L = [l * dPsi / mA for l in range(0, floor(mA * bimp / dPsi) + 1)]
    if abs(L[-1] - bimp) > 1e-12:
        L.append(bimp)
    return L


def block_count(bfun, bimp):
    n = 2  # FA + X2
    for s in "NC":
        n += sum(npts + len(supp(s, J, bfun(J))) for J in SPINS[s])
    return n + 2 + 3 * len(impact_nodes(bimp))  # paRow, pbRow + 3 per impact node


def psi(s, J, M2):
    return (J + 3.5) * theta(s, M2)


# ---------------------------------------------------------------- checks
n20, n3 = block_count(bmax_v20, 100), block_count(bmax_v3, 200)
print("block count test20 grid :", n20, "(1.log reports 61270; 14fa095 script gives this value)")
print("block count test20v3    :", n3, "(1v3.log reports 62301)")
assert n3 == 62301, "v3 recount must match the log"
added = {s: sum(len(supp(s, J, 200)) - len(supp(s, J, 100)) for J in SPINS[s] if J >= 400) for s in "NC"}
print("new supp masses: neutral", added["N"], " charged", added["C"], " (header: 57 x 12 = 684)")
print("impact nodes:", len(impact_nodes(100)), "->", len(impact_nodes(200)),
      " dropped:", sorted(set(impact_nodes(100)) - set(impact_nodes(200)))[-1:])

# phase resolution of the BASE grid in the window 100 < b <= 200 (psi in (100 mA, 200 mA]) for J <= 200
print("\nmax base-grid phase gap inside m_A*100 < psi <= m_A*200  (supp step = pi/4 = 0.785):")
for s in "NC":
    rows = []
    for J in (40, 60, 100, 150, 200):
        if J not in SPINS[s]:
            continue
        ps = sorted(psi(s, J, M2) for M2 in base(s))
        win = [p for p in ps if 100 * mA < p <= 200 * mA]
        gaps = np.diff(sorted(set([p for p in ps if p <= 100 * mA][-1:] + win + [p for p in ps if p > 200 * mA][:1])))
        rows.append(f"J={J}: {len(win)} pts, max gap {gaps.max():.2f}" if len(win) else f"J={J}: window above threshold")
    print(" ", s, "; ".join(rows))

# ---------------------------------------------------------------- plotting helpers
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10})
C_BASE, C_SUPP, C_NEW, C_LINE = "#6b7a8f", "#2a9d8f", "#e63946", "#264653"


def setup(ax, s):
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=200, linscale=1.6)
    ax.set_xlim(0.6, 1e10)
    ax.set_ylim(-5, 8e4)
    ax.set_yticks([0, 50, 100, 150, 200, 400, 1000, 10000, 50000])
    ax.set_yticklabels(["0", "50", "100", "150", "200", "400", "$10^3$", "$10^4$", "$5\\cdot10^4$"])
    ax.axhline(200, color="0.85", lw=0.8, zorder=0)
    ax.axvline(GAP[s], color="k", lw=0.8, ls=":")
    ax.set_xlabel("$M^2$  (units $m_2=1$)")
    ax.set_ylabel("spin $J$")
    ax.set_title(NAME[s])


def scatter_rows(ax, s, fn, **kw):
    xs, ys = [], []
    for J in SPINS[s]:
        m = fn(J)
        xs += list(m)
        ys += [J] * len(m)
    ax.scatter(xs, ys, s=1.2, linewidths=0, **kw)


def b_line(ax, b, **kw):
    M2 = np.logspace(-0.2, 10, 400)
    ax.plot(M2, b * np.sqrt(M2) / 2, **kw)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=170)
    plt.close(fig)
    print("wrote", name)


# Fig 1: base grid only
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, s in zip(axs, "NC"):
    setup(ax, s)
    scatter_rows(ax, s, lambda J: base(s), color=C_BASE, label="base grid (150 Chebyshev masses per spin)")
    ax.legend(loc="lower right", markerscale=6, fontsize=8)
fig.suptitle("Step 1. Base grid: every sampled spin carries the same 150 masses", y=1.0)
save(fig, "fig1_base.png")

# Fig 2: + test20 supplementary masses, b_max = 100
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, s in zip(axs, "NC"):
    setup(ax, s)
    scatter_rows(ax, s, lambda J: base(s), color=C_BASE, label="base grid")
    scatter_rows(ax, s, lambda J: supp(s, J, 100), color=C_SUPP, label="supplementary, $J\\geq40$, $b\\leq100$")
    b_line(ax, 100, color=C_LINE, lw=1, ls="--", label="$b=2J/M=100$")
    ax.axhline(40, color=C_SUPP, lw=0.6, ls=":")
    ax.legend(loc="lower right", markerscale=6, fontsize=8)
fig.suptitle("Step 2. test20.m: phase-spaced supplementary masses fill $b\\lesssim100$ (right of the dashed line)", y=1.0)
save(fig, "fig2_supp_test20.png")

# Fig 3: the v3 additions
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, s in zip(axs, "NC"):
    setup(ax, s)
    scatter_rows(ax, s, lambda J: base(s), color=C_BASE, label="base grid")
    scatter_rows(ax, s, lambda J: supp(s, J, 100), color=C_SUPP, label="supplementary, $b\\leq100$ (test20)")
    scatter_rows(ax, s, lambda J: sorted(set(supp(s, J, bmax_v3(J))) - set(supp(s, J, 100))),
                 color=C_NEW, label="NEW in v3: $J\\geq400$, $100<b\\leq200$")
    b_line(ax, 100, color=C_LINE, lw=1, ls="--", label="$b=100$")
    b_line(ax, 200, color=C_NEW, lw=1, ls="--", label="$b=200$")
    M2 = np.logspace(-0.2, 10, 400)
    ax.fill_betweenx([400, 6e4], 1, 1, alpha=0)  # keep limits
    lo, hi = (2 * 400 / 200) ** 2, 1e10
    ax.axhspan(400, 6e4, xmin=0, xmax=1, color=C_NEW, alpha=0.04, zorder=0)
    ax.axhline(400, color=C_NEW, lw=0.6, ls=":")
    ax.text(2e2, 3000, "A", color=C_NEW, fontsize=12, weight="bold")
    ax.text(4.0, 165, "B", color="k", fontsize=11, weight="bold",
            bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.85))
    ax.legend(loc="lower right", markerscale=6, fontsize=8)
fig.suptitle("Step 3. test20v3.m: for $J\\geq400$ only, supplementary masses continue to $b=200$ (between the dashed lines)", y=1.0)
save(fig, "fig3_v3_new.png")

# Fig 4: phase-coordinate view psi = (J+7/2) theta_s(M) ~ m_A b
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, s in zip(axs, "NC"):
    Js = SPINS[s]
    xb, yb, xs1, ys1, xn, yn = [], [], [], [], [], []
    for J in Js:
        pb = [psi(s, J, M2) for M2 in base(s)]
        xb += [J] * len(pb); yb += pb
        p1 = [psi(s, J, M2) for M2 in supp(s, J, 100)]
        xs1 += [J] * len(p1); ys1 += p1
        pn = [psi(s, J, M2) for M2 in sorted(set(supp(s, J, bmax_v3(J))) - set(supp(s, J, 100)))]
        xn += [J] * len(pn); yn += pn
    ax.scatter(xb, yb, s=1.2, lw=0, color=C_BASE, label="base grid")
    ax.scatter(xs1, ys1, s=1.2, lw=0, color=C_SUPP, label="supplementary $b\\leq100$")
    ax.scatter(xn, yn, s=2.0, lw=0, color=C_NEW, label="NEW in v3 ($J\\geq400$)")
    Jc = np.array(sorted(set(range(0, 201)) | set(Js)), float)
    ax.plot(Jc, [(J + 3.5) * theta(s, GAP[s]) for J in Jc], color="k", lw=0.9, ls=":", label="threshold $M^2=m_{\\rm gap}^2$")
    ax.axhline(100 * mA, color=C_LINE, lw=1, ls="--", label="$\\psi=m_A b$, $b=100$")
    ax.axhline(200 * mA, color=C_NEW, lw=1, ls="--", label="$b=200$")
    ax.set_xscale("symlog", linthresh=200, linscale=1.6)
    ax.set_xticks([0, 50, 100, 150, 200, 400, 1000, 10000, 50000])
    ax.set_xticklabels(["0", "50", "100", "150", "200", "400", "$10^3$", "$10^4$", "$5\\cdot10^4$"])
    ax.set_xlim(-5, 8e4)
    ax.set_ylim(0, 120)
    ax.set_xlabel("spin $J$")
    ax.set_ylabel("$\\psi=(J+7/2)\\,\\theta_s(M)\\;\\approx\\; m_A b$")
    ax.set_title(NAME[s])
    ax.set_ylim(0, 135)
    ax.legend(loc="upper right", markerscale=6, fontsize=7.5, framealpha=0.95)
    ax.text(1500, 67, "A", color=C_NEW, fontsize=12, ha="center", weight="bold")
    for xl, yl, lab in ((150, 67, "B"), (150, 22, "C")):
        ax.text(xl, yl, lab, color="k", fontsize=11, ha="center", weight="bold",
                bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.85))
fig.suptitle("Step 4. Phase coordinate $\\psi\\approx m_A b$: supplementary rows at $\\psi=j\\pi/4$; "
             "v3 lifts the ceiling to $m_A\\cdot200$ for $J\\geq400$ only", y=1.0)
save(fig, "fig4_phase_view.png")

# Fig 5: impact-parameter nodes
fig, ax = plt.subplots(figsize=(11, 3.4))
n20b, n3b = impact_nodes(100), impact_nodes(200)
ax.scatter(n20b[:-1], [1] * (len(n20b) - 1), s=14, color=C_LINE, label=f"test20: {len(n20b)} nodes, $0\leq b\leq100$")
ax.scatter([100], [1], s=70, facecolors="none", edgecolors=C_LINE, lw=1.5,
           label="$b=100$ endpoint node: test20 only (nearest v3 nodes: %.2f, %.2f)" % (56 * dPsi / mA, 57 * dPsi / mA))
keep = [b for b in n3b if b <= 100]
new = [b for b in n3b if b > 100]
ax.scatter(keep, [0] * len(keep), s=14, color=C_SUPP, label="test20v3: nodes shared with test20")
ax.scatter(new, [0] * len(new), s=14, color=C_NEW, label=f"test20v3: {len(new)} nodes with $100<b\leq200$ (incl. endpoint 200)")
ax.axvline(100, color="0.6", lw=0.8, ls="--")
ax.set_yticks([0, 1]); ax.set_yticklabels(["test20v3", "test20"])
ax.set_ylim(-0.6, 1.6); ax.set_xlim(-3, 205)
ax.set_xlabel("impact parameter $b=2J/M$   (node spacing $\Delta\psi/m_A=%.3f$)" % (dPsi / mA))
ax.legend(loc="upper center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, -0.32), frameon=False)
ax.set_title(r"Step 5. Fixed-impact-parameter ($J,M\to\infty$, $b$ fixed) blocks: 3 blocks per node")
save(fig, "fig5_impact_nodes.png")
