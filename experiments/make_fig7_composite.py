#!/usr/bin/env python3
"""
Composite Fig. 7 (Independent-engine validation, E6) for the paper.

Assembles two sub-panels into one figure, matplotlib-styled to match the
rest of the paper's figures (make_figures.py):
  (a) Simscape Multibody stick/slide sweep at the central (nominal) realization
  (b) Simscape mu-sweep recovering atan(mu)

A PyBullet panel is not included: its "predicted=observed" residual is not
used as validation evidence (see
supplementary/pybullet_identification_retraction.md). PyBullet's raw,
unfitted transition (E6(ii)) and the box-spanning PyBullet sweep (E6(iii))
are reported only in text, not in this figure.

This script was reconstructed because the original composite-assembly script
was not present in the reproducibility package (only the two underlying
per-engine sweep scripts were: supplementary/robust_motion_cone_simscape.m,
supplementary/cone_simscape_sweep.m).
Data sources:
  (a), (b) -- parsed from the console logs of the two .m scripts (they do not
             export a machine-readable results file; see the parsing helpers
             below). Run robust_motion_cone_simscape.m and
             cone_simscape_sweep.m in MATLAB first and pass the log paths.

Usage:
    python make_fig7_composite.py --main-log <path> --mu-log <path>
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT_FIGS = ROOT / "figures"

PUSHER_SPEED = 0.1
VREL_TOL = 0.15 * PUSHER_SPEED  # matches robust_motion_cone_simscape.m


def load_main_sweep_json(path: Path):
    """Preferred data source: supplementary/main_sweep_results.json, written
    directly by robust_motion_cone_simscape.m (no log-parsing needed)."""
    d = json.loads(path.read_text())
    return {
        "angles": np.array(d["angle_deg"]), "sticking": np.array(d["sticking"], dtype=bool),
        "nom_lo": d["nom_lo_deg"], "nom_hi": d["nom_hi_deg"],
        "rob_lo": d["rob_lo_deg"], "rob_hi": d["rob_hi_deg"],
        "accuracy": d["accuracy_pct"], "n_correct": d["n_correct"], "n_total": d["n_total"],
        "precision": d["precision_pct"], "recall": d["recall_pct"], "mu0": d["mu0"],
    }


def load_mu_sweep_json(path: Path):
    """Preferred data source: supplementary/mu_sweep_results.json, written
    directly by cone_simscape_sweep.m (no log-parsing needed)."""
    d = json.loads(path.read_text())
    return {"mu": np.array(d["mu"]), "theory": np.array(d["theory_deg"]),
            "sim": np.array(d["sim_deg"]), "err": np.array(d["err_deg"])}


def parse_main_sweep_log(path: Path):
    text = path.read_text(errors="replace")
    angles, vrel_t = [], []
    for m in re.finditer(
        r"^\s*([+-]?\d+) deg \| Vx=[-\d.]+ Vy=[-\d.]+ dir=[-\d.]+ "
        r"vy_rat=[-\d.]+ vrel_t=([-\d.]+)", text, re.MULTILINE):
        angles.append(int(m.group(1)))
        vrel_t.append(float(m.group(2)))
    if not angles:
        raise ValueError(f"No per-angle rows parsed from {path}")

    def _f(pat, cast=float):
        mo = re.search(pat, text)
        if not mo:
            raise ValueError(f"Pattern not found in {path}: {pat}")
        return cast(mo.group(1))

    nom_lo = _f(r"Nominal cone:\s*\[\s*([-\d.]+),")
    nom_hi = _f(r"Nominal cone:\s*\[\s*[-\d.]+,\s*([-\d.]+)\]")
    rob_lo = _f(r"Robust cone:\s*\[\s*([-\d.]+),")
    rob_hi = _f(r"Robust cone:\s*\[\s*[-\d.]+,\s*([-\d.]+)\]")
    accuracy = _f(r"Overall accuracy:\s*([\d.]+)%")
    acc_frac = re.search(r"Overall accuracy:\s*[\d.]+%\s*\((\d+)/(\d+)\)", text)
    precision = _f(r"Precision \(stick\):\s*([\d.]+)%")
    recall = _f(r"Recall \(stick\):\s*([\d.]+)%")
    mu0 = _f(r"Contact friction: mu = ([\d.]+)")

    angles = np.array(angles)
    vrel_t = np.array(vrel_t)
    sticking = np.abs(vrel_t) < VREL_TOL
    return {
        "angles": angles, "sticking": sticking,
        "nom_lo": nom_lo, "nom_hi": nom_hi, "rob_lo": rob_lo, "rob_hi": rob_hi,
        "accuracy": accuracy, "n_correct": int(acc_frac.group(1)),
        "n_total": int(acc_frac.group(2)),
        "precision": precision, "recall": recall, "mu0": mu0,
    }


def parse_mu_sweep_log(path: Path):
    text = path.read_text(errors="replace")
    rows = re.findall(
        r"^\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s*$", text, re.MULTILINE)
    rows = [r for r in rows if 0 < float(r[0]) < 1]  # keep only mu rows (0<mu<1)
    if not rows:
        raise ValueError(f"No mu-sweep summary rows parsed from {path}")
    mu = np.array([float(r[0]) for r in rows])
    theory = np.array([float(r[1]) for r in rows])
    sim = np.array([float(r[2]) for r in rows])
    err = np.array([float(r[3]) for r in rows])
    return {"mu": mu, "theory": theory, "sim": sim, "err": err}


plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 11,
    "axes.linewidth": 0.6, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def _stick_slide_panel(ax, angles, sticking, edges, edge_style, title):
    """One (a)/(c)-style row: green dots = stick, red squares = slide,
    shaded band + solid lines at `edges`, dashed/dotted reference lines from
    `edge_style` = list of (value, color, linestyle)."""
    stick_a = angles[sticking]
    slide_a = angles[~sticking]
    ax.scatter(stick_a, np.ones_like(stick_a), marker="o", s=45,
               color="#2a8c4a", zorder=3)
    ax.scatter(slide_a, np.zeros_like(slide_a), marker="s", s=55,
               color="#c0392b", zorder=3)
    lo, hi = edges
    ax.axvspan(lo, hi, color="#2a8c4a", alpha=0.12, zorder=1)
    ax.axvline(lo, color="#1a1a1a", lw=1.6, zorder=2)
    ax.axvline(hi, color="#1a1a1a", lw=1.6, zorder=2)
    for v, c, ls in edge_style:
        ax.axvline(v, color=c, lw=1.1, ls=ls, zorder=2)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["slide", "stick"])
    ax.set_ylim(-0.5, 1.5)
    ax.set_xlim(-90, 90)
    ax.set_xlabel("push angle (deg)")
    ax.set_title(title, fontsize=10)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)


def make_composite(main_src: Path, mu_src: Path, out_name="fig5_engine_sweep_validated"):
    """main_src / mu_src may be either the .m scripts' JSON exports
    (supplementary/main_sweep_results.json, mu_sweep_results.json -- preferred)
    or their console-log captures (fallback, used if JSON wasn't exported).

    NOTE: this composite renders 2 panels (Simscape stick/slide sweep +
    Simscape mu-sweep). A PyBullet panel is not included: its
    "predicted=observed" residual is not used as validation evidence (see
    supplementary/pybullet_identification_retraction.md); PyBullet's raw,
    unfitted transition is reported only in text (E6(ii)), not in this
    figure. pybullet_sweep.py still produces its own standalone figure when
    run; it is not part of this composite and is not tracked in the repo."""
    a = load_main_sweep_json(main_src) if main_src.suffix == ".json" else parse_main_sweep_log(main_src)
    b = load_mu_sweep_json(mu_src) if mu_src.suffix == ".json" else parse_mu_sweep_log(mu_src)

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.6))

    # (a) Simscape Multibody stick/slide sweep
    _stick_slide_panel(
        axes[0], a["angles"], a["sticking"], (a["rob_lo"], a["rob_hi"]),
        [(a["nom_lo"], "#7a7a7a", "--"), (a["nom_hi"], "#7a7a7a", "--")],
        f"(a) Simscape Multibody stick/slide sweep  (central realization "
        f"$\\mu_0$={a['mu0']:.2f}, $c_0$=0.06; {a['accuracy']:.1f}% accuracy, "
        f"{a['n_correct']}/{a['n_total']})",
    )

    # (b) mu-sweep
    ax = axes[1]
    mu_fine = np.linspace(min(b["mu"]) - 0.03, max(b["mu"]) + 0.03, 200)
    ax.plot(mu_fine, np.rad2deg(np.arctan(mu_fine)), "k--", lw=1.6, label=r"atan($\mu$)")
    ax.plot(b["mu"], b["sim"], "o", color="#2a6fc0", ms=9, mfc="none", mew=1.8,
            label="Simscape")
    ax.set_xlabel(r"friction coefficient $\mu$")
    ax.set_ylabel("cone half-angle (deg)")
    max_err = float(np.max(np.abs(b["err"])))
    ax.set_title(f"(b)  Simscape $\\mu$-sweep recovers atan($\\mu$) exactly  "
                 f"(max err $\\leq$ {max_err:.2f}$^\\circ$)", fontsize=10)
    ax.legend(loc="upper left", fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.tight_layout(h_pad=2.2)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_FIGS / f"{out_name}.{ext}")
    plt.close(fig)
    print(f"Wrote {OUT_FIGS / out_name}.png/.pdf")
    print(f"  (a) accuracy={a['accuracy']}% precision={a['precision']}% "
          f"recall={a['recall']}% robust=[{a['rob_lo']},{a['rob_hi']}] "
          f"nominal=[{a['nom_lo']},{a['nom_hi']}]")
    print(f"  (b) max_err={max_err:.4f} deg  mu={list(b['mu'])}")
    print("  (Note: PyBullet panel (c) is not part of this composite; see "
          "supplementary/pybullet_identification_retraction.md. "
          "pybullet_sweep.py produces its own standalone figure when run.)")


if __name__ == "__main__":
    SUPP = ROOT / "supplementary"
    ap = argparse.ArgumentParser()
    ap.add_argument("--main-log", type=Path, default=SUPP / "main_sweep_results.json",
                     help="main_sweep_results.json (preferred) or a console-log capture")
    ap.add_argument("--mu-log", type=Path, default=SUPP / "mu_sweep_results.json",
                     help="mu_sweep_results.json (preferred) or a console-log capture")
    args = ap.parse_args()
    make_composite(args.main_log, args.mu_log)
