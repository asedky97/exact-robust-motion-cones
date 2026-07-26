#!/usr/bin/env python3
"""
Box-spanning PyBullet validation: tests the robust cone itself against an
independent solver across the full uncertainty box, rather than a single
fitted point estimate.

Unlike pybullet_sweep.py (E6), which identifies ONE effective (mu_eff, c_eff)
realization and checks predicted==observed at that single point, this script
tests the ROBUST CONE ITSELF against an independent solver by:
  1. Sampling mu uniformly across the FULL uncertainty box [0.25, 0.65]
     (directly controllable in PyBullet via changeDynamics), not a single
     fitted value.
  2. For each realized mu, commanding a body-frame push direction sampled
     from three zones -- inside the robust cone [Phi-, Phi+], inside nominal
     but outside robust, and outside nominal -- mirroring the paper's E2
     zone protocol, but in an independent rigid-body engine instead of the
     analytical resolve_contact model.
  3. Verifying zero slides for Zone R (the actual claim under test: does the
     analytically-worst-case-safe cone hold up against an independent
     contact solver as mu ranges over the box, not just at one point).

Limitation (stated explicitly, not hidden): PyBullet's rigid-body primitive
does not expose an independently controllable "c" (limit-surface length) or
"alpha" (contact-normal angle) the way it exposes lateralFriction for mu, so
this test spans the box in mu only; c is fixed at the block's true rigid-body
value (as in pybullet_sweep.py's single-realization E6 test). This is the
same c-limitation already stated for the Simscape/PyBullet single-point test
in supplementary/pybullet_identification_retraction.md.

The board's half-extents (pybullet_box_sweep_world.py's BOARD_HALF) are
chosen so its implied c -- the radius of gyration of its footprint,
sqrt((a^2+b^2)/3), the same formula the paper uses for a uniform rectangular
support (Sec. II.A) -- falls inside the assumed c-box [0.045, 0.075] m,
matching the model this test is checking. The contact point (CONTACT_X in
pybullet_box_sweep_world.py) sits at the board's edge, x = -a, consistent
with ContactGeometry(p0=(CONTACT_X, 0)) below.

Usage:
    python pybullet_box_sweep.py [--n 300]
Output:
    experiments/pybullet_box_sweep_results.json
    ../figures/fig_pybullet_box_sweep.png/.pdf
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import ContactGeometry, UncertaintyBox, cone_angles, robust_cone_exact  # noqa

from pybullet_box_sweep_world import build_world, run_push_mu, CONTACT_X  # noqa


def _wilson_ci(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return (float(round(max(0.0, 100 * (centre - margin)), 2)),
            float(round(min(100.0, 100 * (centre + margin)), 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300, help="trials per zone")
    ap.add_argument("--gui", action="store_true")
    args = ap.parse_args()

    box = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                         alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))
    geom = ContactGeometry(p0=np.array([CONTACT_X, 0.0]), n0=np.array([1.0, 0.0]))
    lo, hi, ref, _ = robust_cone_exact(geom, box)
    mu0, c0, a0, d0 = box.nominal()
    lon, hin = cone_angles(geom, mu0, c0, a0, d0, ref)
    ROB_LO, ROB_HI = np.rad2deg(lo), np.rad2deg(hi)
    NOM_LO, NOM_HI = np.rad2deg(lon), np.rad2deg(hin)
    print(f"Robust cone [{ROB_LO:.2f}, {ROB_HI:.2f}]  Nominal [{NOM_LO:.2f}, {NOM_HI:.2f}]")

    cid, board, pusher = build_world(args.gui)
    rng = np.random.default_rng(42)

    zones = {"R": {"lo": ROB_LO, "hi": ROB_HI, "slides": 0, "n": 0},
             "N": {"lo": None, "hi": None, "slides": 0, "n": 0},  # nominal-only band, sampled explicitly below
             "S": {"lo": None, "hi": None, "slides": 0, "n": 0}}

    t0 = time.perf_counter()
    records = []
    for zone in ("R", "N", "S"):
        for _ in range(args.n):
            mu = rng.uniform(*box.mu)  # span the FULL box, independent-solver test
            if zone == "R":
                psi = rng.uniform(ROB_LO, ROB_HI)
            elif zone == "N":
                # inside nominal, outside robust (two side-bands)
                side = rng.choice([-1, 1])
                if side < 0:
                    psi = rng.uniform(NOM_LO, ROB_LO)
                else:
                    psi = rng.uniform(ROB_HI, NOM_HI)
            else:  # S: outside nominal
                side = rng.choice([-1, 1])
                if side < 0:
                    psi = rng.uniform(NOM_LO - 20, NOM_LO)
                else:
                    psi = rng.uniform(NOM_HI, NOM_HI + 20)

            mode, eff = run_push_mu(board, pusher, psi, mu)
            slide = mode != "stick"
            zones[zone]["slides"] += int(slide)
            zones[zone]["n"] += 1
            records.append({"zone": zone, "mu": float(mu), "psi_deg": float(psi),
                            "mode": mode, "efficiency": float(eff)})
    elapsed = time.perf_counter() - t0
    print(f"Swept {len(records)} trials in {elapsed:.1f}s")

    res = {"robust_deg": [round(ROB_LO, 2), round(ROB_HI, 2)],
           "nominal_deg": [round(NOM_LO, 2), round(NOM_HI, 2)],
           "mu_box": list(box.mu)}
    for z in ("R", "N", "S"):
        k, n = zones[z]["slides"], zones[z]["n"]
        ci = _wilson_ci(k, n)
        res[z] = {"slide_pct": round(100 * k / max(n, 1), 2), "n": n, "ci_95": ci}
        print(f"  Zone {z}: {res[z]['slide_pct']}% slide  (n={n}, 95% CI {ci})")

    (ROOT / "experiments" / "pybullet_box_sweep_results.json").write_text(
        json.dumps({"summary": res, "records": records}, indent=2))
    print(f"Wrote {ROOT / 'experiments' / 'pybullet_box_sweep_results.json'}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "cm", "font.size": 11,
        "axes.linewidth": 0.6, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    fig, ax = plt.subplots(figsize=(4, 3))
    zones_lbl = ["Zone R\n(robust)", "Zone N\n(nominal)", "Zone S\n(outside)"]
    vals = [res["R"]["slide_pct"], res["N"]["slide_pct"], res["S"]["slide_pct"]]
    cols = ["#2a8c4a", "#cc8800", "#c0392b"]
    bars = ax.bar(zones_lbl, vals, color=cols, width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, min(v, 95) + 2, f"{v:.2f}%",
                ha="center", fontsize=9, fontweight="bold")
    ax.set_ylabel("sliding failure rate (%)"); ax.set_ylim(0, 110)
    ax.set_title(f"PyBullet, $\\mu$ spanning box (n={args.n}/zone)", fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "fig_pybullet_box_sweep.png")
    fig.savefig(ROOT / "figures" / "fig_pybullet_box_sweep.pdf")
    print(f"Wrote {ROOT / 'figures' / 'fig_pybullet_box_sweep.png'}")


if __name__ == "__main__":
    main()
