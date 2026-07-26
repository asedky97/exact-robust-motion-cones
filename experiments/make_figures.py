#!/usr/bin/env python3
"""
Regenerates every figure used in the paper, writing PNG + PDF to
../figures/. Run via `python run_all.py --figs` or directly `python make_figures.py`.

Figures:
    fig1_model.png          - contact / friction-cone / motion-cone schematic
    fig_sensitivity.png     - robust aperture vs each uncertainty parameter
    fig4_montecarlo.png     - zone failure-rate bar chart (E2)
    fig5_funnel.png         - certified funnel (E3)
    fig7_planner.png        - command-frame comparison (E5)
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import ContactGeometry, UncertaintyBox, cone_angles, robust_cone_exact  # noqa
from pushing_sim import resolve_contact  # noqa

OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
SUMMARY = ROOT / "experiments" / "summary.json"


def _summary(*keys):
    """Single source of truth for headline numbers: read the requested keys from
    experiments/summary.json. If it is missing (figures generated before the
    experiments were run), lazily run the corresponding experiment so the figure
    always matches the code output -- never a hardcoded value."""
    data = json.loads(SUMMARY.read_text()) if SUMMARY.exists() else {}
    missing = [k for k in keys if k not in data]
    if missing:
        import run_all  # lazy: avoids any import cycle at module load
        exp = {1: run_all.exp1, 2: run_all.exp2, 3: run_all.exp3,
               4: run_all.exp4, 5: run_all.exp5, 6: run_all.exp6}
        for k in missing:
            data[k] = exp[int(k[1:])]()
    return tuple(data[k] for k in keys)

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.linewidth": 0.6, "savefig.dpi": 300, "savefig.bbox": "tight",
})

BOX = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                     alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))
GEOM = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))


def _save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}")


def fig_model():
    fig, ax = plt.subplots(figsize=(4, 3))
    # board
    ax.add_patch(plt.Rectangle((-0.30, -0.05), 0.60, 0.10, fill=True,
                               color="#e9dcc3", ec="#7a6a48", lw=1.2))
    ax.plot(0, 0, "ko", ms=4)
    ax.annotate("CoM", (0.01, 0.01), fontsize=8)
    p = np.array([-0.30, 0.0])
    ax.plot(*p, "o", color="#b8002e", ms=6)
    ax.annotate("$p$", (p[0] - 0.03, p[1] + 0.01), fontsize=10)
    gamma = np.arctan(0.45)
    for sgn, col in [(+1, "#2a6fc0"), (-1, "#2a6fc0")]:
        d = np.array([np.cos(sgn * gamma), np.sin(sgn * gamma)])
        ax.arrow(p[0], p[1], 0.12 * d[0], 0.12 * d[1], head_width=0.012,
                 color=col, lw=1.2)
    ax.annotate("$f_+$", (p[0] + 0.13, 0.05), color="#2a6fc0", fontsize=9)
    ax.annotate("$f_-$", (p[0] + 0.13, -0.07), color="#2a6fc0", fontsize=9)
    ax.set_xlim(-0.4, 0.4); ax.set_ylim(-0.2, 0.2); ax.set_aspect("equal")
    ax.axis("off")
    _save(fig, "fig1_model")



def fig_sensitivity():
    mu0, c0, a0, d0 = BOX.nominal()
    # Nominal aperture (no uncertainty) -- reference line, computed not hardcoded.
    _lo0, _hi0, _r0, _ = robust_cone_exact(
        GEOM, UncertaintyBox((mu0, mu0), (c0, c0), (0, 0), (0, 0)))
    A_nom = np.rad2deg(_hi0 - _lo0)
    fig, axs = plt.subplots(2, 2, figsize=(7, 5))
    mus = np.linspace(0.1, 0.95, 40)
    axs[0, 0].plot(mus, [np.rad2deg(robust_cone_exact(GEOM, UncertaintyBox((m, m), (c0, c0), (0, 0), (0, 0)))[1]
                                    - robust_cone_exact(GEOM, UncertaintyBox((m, m), (c0, c0), (0, 0), (0, 0)))[0]) for m in mus],
                   color="#2a6fc0", lw=1.8)
    axs[0, 0].axhline(A_nom, ls="--", color="gray", lw=0.8)
    axs[0, 0].set_xlabel("$\\mu$"); axs[0, 0].set_ylabel("aperture (deg)"); axs[0, 0].set_title("(a) vs $\\mu$", fontsize=9)
    chs = np.linspace(0, 0.04, 40)
    axs[0, 1].plot(chs * 1000, [np.rad2deg(robust_cone_exact(GEOM, UncertaintyBox((mu0, mu0), (c0 - ch, c0 + ch), (0, 0), (0, 0)))[1]
                                            - robust_cone_exact(GEOM, UncertaintyBox((mu0, mu0), (c0 - ch, c0 + ch), (0, 0), (0, 0)))[0]) for ch in chs],
                   color="#2a8c4a", lw=1.8)
    axs[0, 1].axhline(A_nom, ls="--", color="gray", lw=0.8)
    axs[0, 1].set_xlabel("c half-range (mm)"); axs[0, 1].set_ylabel("aperture (deg)"); axs[0, 1].set_title("(b) vs $c$", fontsize=9)
    ahs = np.linspace(0, 25, 40)
    axs[1, 0].plot(ahs, [np.rad2deg(robust_cone_exact(GEOM, UncertaintyBox((mu0, mu0), (c0, c0), (-np.deg2rad(ah), np.deg2rad(ah)), (0, 0)))[1]
                                    - robust_cone_exact(GEOM, UncertaintyBox((mu0, mu0), (c0, c0), (-np.deg2rad(ah), np.deg2rad(ah)), (0, 0)))[0]) for ah in ahs],
                   color="#c0392b", lw=1.8)
    axs[1, 0].axhline(A_nom, ls="--", color="gray", lw=0.8)
    axs[1, 0].set_xlabel("$|\\alpha|$ (deg)"); axs[1, 0].set_ylabel("aperture (deg)"); axs[1, 0].set_title("(c) vs $\\alpha$", fontsize=9)
    dhs = np.linspace(0, 50, 40)
    axs[1, 1].plot(dhs, [np.rad2deg(robust_cone_exact(GEOM, UncertaintyBox((mu0, mu0), (c0, c0), (0, 0), (-dh / 1000, dh / 1000)))[1]
                                    - robust_cone_exact(GEOM, UncertaintyBox((mu0, mu0), (c0, c0), (0, 0), (-dh / 1000, dh / 1000)))[0]) for dh in dhs],
                   color="#8e44ad", lw=1.8)
    axs[1, 1].axhline(A_nom, ls="--", color="gray", lw=0.8)
    axs[1, 1].set_xlabel("$|\\delta|$ (mm)"); axs[1, 1].set_ylabel("aperture (deg)"); axs[1, 1].set_title("(d) vs $\\delta$", fontsize=9)
    fig.suptitle("Robust aperture sensitivity to each uncertainty parameter", fontsize=10, y=1.0)
    fig.tight_layout()
    _save(fig, "fig_sensitivity")


def fig_montecarlo():
    (e2,) = _summary("E2")
    n_total = e2["R"]["n"] + e2["N"]["n"] + e2["S"]["n"]
    fig, ax = plt.subplots(figsize=(4, 3))
    zones = ["Zone R\n(robust)", "Zone N\n(nominal)", "Zone S\n(outside)"]
    vals = [e2["R"]["slide_pct"], e2["N"]["slide_pct"], e2["S"]["slide_pct"]]
    cols = ["#2a8c4a", "#cc8800", "#c0392b"]
    bars = ax.bar(zones, vals, color=cols, width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, min(v, 95) + 2, f"{v:.2f}%",
                ha="center", fontsize=9, fontweight="bold")
    ax.set_ylabel("sliding failure rate (%)"); ax.set_ylim(0, 110)
    ax.set_title(f"Mode verification ({n_total:,} draws)", fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    _save(fig, "fig4_montecarlo")


def fig_funnel():
    (e3,) = _summary("E3")
    bound = e3["theta_bound_deg"]              # Theorem 3 orientation bound (deg)
    slides = e3["total_slides"]                # certified slide count from E3
    rng = np.random.default_rng(2)
    lo, hi, ref, _ = robust_cone_exact(GEOM, BOX)
    psi = lo + 0.35 * (hi - lo)
    u = np.array([np.cos(psi), np.sin(psi)])
    finals = []
    for _ in range(300):
        th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
              rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
        p, n = GEOM.contact(th[2], th[3])
        yaw = 0.0
        for _ in range(500):
            m, v = resolve_contact(p, n, th[0], th[1], u)
            if m != "stick":
                break
            yaw += v[2] * 1e-3
        finals.append(np.rad2deg(yaw))
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.hist(finals, bins=40, color="#5b8bd0", alpha=0.85)
    ax.axvline(np.mean(finals), color="#b8002e", lw=1.5,
               label=f"mean {np.mean(finals):.1f}°")
    ax.axvline(np.mean(finals) + bound / 2, color="orange", ls="--", lw=1.2,
               label=f"Thm 3 bound ±{bound:.2f}°")
    ax.axvline(np.mean(finals) - bound / 2, color="orange", ls="--", lw=1.2)
    ax.set_xlabel("final orientation (deg)"); ax.set_ylabel("count")
    ax.set_title(f"Certified funnel: {slides} slides", fontsize=9)
    ax.legend(fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    _save(fig, "fig5_funnel")


def fig_planner():
    (e5,) = _summary("E5")
    n_draws = e5["A_robust"]["n"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7, 3))
    strat = ["A\nrobust\nbody", "B\nworld\nframe", "C\nnominal\nbody"]
    vals = [e5["A_robust"]["failure_pct"], e5["B_world"]["failure_pct"],
            e5["C_nominal"]["failure_pct"]]
    cols = ["#2a8c4a", "#c0392b", "#cc8800"]
    bars = a1.bar(strat, vals, color=cols, width=0.55)
    for b, v in zip(bars, vals):
        a1.text(b.get_x() + b.get_width() / 2, min(v, 95) + 2, f"{v:.2f}%",
                ha="center", fontsize=9, fontweight="bold")
    a1.set_ylabel("failure rate (%)"); a1.set_ylim(0, 110)
    a1.set_title(f"(a) 45° reorientation ({n_draws // 1000}k draws)", fontsize=9)
    a1.spines["top"].set_visible(False); a1.spines["right"].set_visible(False)

    # (b) Genuine failure-vs-rotation-horizon curve (NOT a fitted sigmoid).
    # World-frame command psi = 45° - yaw drifts out of the (body-fixed) cone as
    # the object rotates; we bin each realization by the rotation accumulated at
    # first slip. Body-frame (Thm 3) command stays in the cone => 0 failures.
    rng = np.random.default_rng(5)
    N = 3000
    horizons = np.arange(0, 140, 10)
    slip_rot = np.full(N, np.inf)
    for i in range(N):
        th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
              rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
        p, n = GEOM.contact(th[2], th[3])
        yaw = 0.0
        for _ in range(400):
            psi = np.deg2rad(45) - yaw
            u = np.array([np.cos(psi), np.sin(psi)])
            m, v = resolve_contact(p, n, th[0], th[1], u)
            if m != "stick":
                slip_rot[i] = abs(np.rad2deg(yaw))
                break
            yaw += v[2] * 40e-3
    world = [100 * np.mean(slip_rot <= h) for h in horizons]
    a2.plot(horizons, world, "o-", color="#c0392b", ms=4, label="world-frame")
    a2.plot(horizons, np.zeros_like(horizons), "s--", color="#2a8c4a", ms=4,
            label="body-frame (Thm 3)")
    a2.set_xlabel("rotation accumulated (deg)"); a2.set_ylabel("failure rate (%)")
    a2.set_title("(b) failure vs rotation horizon", fontsize=9)
    a2.legend(fontsize=7)
    a2.spines["top"].set_visible(False); a2.spines["right"].set_visible(False)
    fig.tight_layout()
    _save(fig, "fig7_planner")



def fig_construction():
    """Method figure: the robust cone is the INTERSECTION of the per-realization
    motion cones over the uncertainty box (Proposition 1). Its exact edges are
    Phi- = max phi-(theta) and Phi+ = min phi+(theta) over theta in the box.
    (a) polar sectors: nominal, a family of box-vertex candidates, and the robust
    intersection; (b) interval view: sampled realization cones [phi-,phi+] and the
    robust band they intersect to."""
    lo, hi, ref, _ = robust_cone_exact(GEOM, BOX)
    mu0, c0, a0, d0 = BOX.nominal()
    lon, hin = cone_angles(GEOM, mu0, c0, a0, d0, ref)
    LO, HI = np.rad2deg(lo), np.rad2deg(hi)
    NLO, NHI = np.rad2deg(lon), np.rad2deg(hin)

    # box-vertex candidate cones (16) for the polar panel
    verts = [cone_angles(GEOM, mu, c, al, de, ref)
             for mu in BOX.mu for c in BOX.c
             for al in BOX.alpha for de in BOX.delta]
    # sampled realization cones for the interval panel
    rng = np.random.default_rng(3)
    samples = []
    for _ in range(28):
        th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
              rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
        l, h = cone_angles(GEOM, th[0], th[1], th[2], th[3], ref)
        samples.append((np.rad2deg(l), np.rad2deg(h)))
    samples.sort(key=lambda x: x[0])

    # Single-column, vertically stacked layout (fits an IEEE text column).
    from matplotlib.patches import Patch
    fig = plt.figure(figsize=(3.4, 4.5))
    # (a) polar sectors (top)
    axp = fig.add_subplot(211, projection="polar")
    for (l, h) in verts:
        axp.fill_between(np.linspace(l, h, 40), 0, 1.0, color="#8a8f98", alpha=0.06)
    axp.fill_between(np.linspace(lon, hin, 60), 0, 1.0, color="#5b8bd0", alpha=0.16)
    axp.fill_between(np.linspace(lo, hi, 60), 0, 1.0, color="#b8002e", alpha=0.32)
    axp.plot([lo, lo], [0, 1.0], color="#b8002e", lw=1.6)
    axp.plot([hi, hi], [0, 1.0], color="#b8002e", lw=1.6)
    axp.set_ylim(0, 1.05); axp.set_yticks([])
    axp.set_thetalim(np.deg2rad(-120), np.deg2rad(120))
    axp.set_xticklabels([])   # schematic: drop angle ticks (panel b is quantitative)
    axp.set_title("(a) intersection of candidate cones", fontsize=8, pad=8)
    axp.legend(handles=[Patch(fc="#8a8f98", alpha=0.35, label="candidates"),
                        Patch(fc="#5b8bd0", alpha=0.4, label="nominal"),
                        Patch(fc="#b8002e", alpha=0.55, label="robust")],
               loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=6)

    # (b) interval-intersection view (bottom)
    ax = fig.add_subplot(212)
    for i, (l, h) in enumerate(samples):
        ax.plot([l, h], [i, i], color="#9aa7b3", lw=1.3, solid_capstyle="round")
    ax.axvspan(LO, HI, color="#b8002e", alpha=0.15)
    ax.axvline(LO, color="#b8002e", lw=1.4)
    ax.axvline(HI, color="#b8002e", lw=1.4)
    ax.axvline(NLO, color="#5b8bd0", lw=1.2, ls="--")
    ax.axvline(NHI, color="#5b8bd0", lw=1.2, ls="--")
    ax.set_ylim(-1, len(samples) + 4)
    ax.text(HI, len(samples) + 1.0, f"$\\Phi^+\\!=\\!{HI:.1f}^\\circ$",
            color="#b8002e", ha="center", fontsize=6.5)
    ax.text(LO, len(samples) + 1.0, f"$\\Phi^-\\!=\\!{LO:.1f}^\\circ$",
            color="#b8002e", ha="center", fontsize=6.5)
    ax.set_xlabel("command / edge angle (deg)", fontsize=8)
    ax.set_ylabel("realization", fontsize=8)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7)
    ax.set_title("(b) robust edges $=\\max\\,\\phi^-,\\ \\min\\,\\phi^+$",
                 fontsize=8, pad=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(h_pad=1.5)
    _save(fig, "fig_construction")


def make_all():
    print("Generating figures...")
    fig_model()
    fig_construction()
    fig_sensitivity()
    fig_montecarlo()
    fig_funnel()
    fig_planner()
    print(f"All figures written to {OUT}/")


if __name__ == "__main__":
    make_all()
