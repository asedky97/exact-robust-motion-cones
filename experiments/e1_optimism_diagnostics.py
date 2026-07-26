#!/usr/bin/env python3
"""
Reproduces three E1 sub-claims (Section V.A, baseline-comparison paragraph)
that are not covered by run_all.py's exp1, which only reports the two
headline numbers (pct_optimistic, max_vtx_err_deg):

  1. The median vertex-only optimism error ("median error 0.04 deg").
  2. That optimistic cases correlate with longer lever arms and more
     off-axis geometry than the sampled population as a whole.
  3. The delta-box-width sensitivity ("a wider delta-box grows both
     incidence and magnitude").

Uses the exact same geometry sampling and vertex-only-vs-exact comparison
as run_all.py::exp1 (verified to reproduce its 1.85% / 0.24 deg headline
numbers at --mult 1.0).

Usage:
    python e1_optimism_diagnostics.py [--mult 3.0] [--n 2000]
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import (  # noqa
    ContactGeometry, UncertaintyBox, robust_cone_exact, cone_angles,
)

DELTA0 = 0.015  # paper-scenario delta half-width [m]


def make_box(delta_mult):
    d = DELTA0 * delta_mult
    return UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                          alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-d, d))


def optimism_sweep(box, n_geom, seed=0):
    """Same geometry sampling and vertex-only-vs-exact comparison as
    run_all.py::exp1's optimism check, additionally recording each sampled
    geometry's lever arm and |off-axis angle| (for the correlation claim)."""
    rng = np.random.default_rng(seed)
    all_r, all_ang = [], []
    opt_r, opt_ang, opt_err = [], [], []
    for _ in range(n_geom):
        r = rng.uniform(0.20, 0.35)
        th_p = rng.uniform(-0.15, 0.15)
        p = np.array([-r * np.cos(th_p), r * np.sin(th_p)])
        ang = rng.uniform(-0.15, 0.15)
        all_r.append(r)
        all_ang.append(abs(ang))
        G = ContactGeometry(p0=p, n0=np.array([np.cos(ang), np.sin(ang)]))
        le, he, r2, _ = robust_cone_exact(G, box)
        if he <= le:
            continue
        vlo, vhi = -np.inf, np.inf
        for mu in box.mu:
            for c in box.c:
                for al in box.alpha:
                    for de in box.delta:
                        l, h = cone_angles(G, mu, c, al, de, r2)
                        vlo = max(vlo, l)
                        vhi = min(vhi, h)
        if vlo < le - 1e-9 or vhi > he + 1e-9:
            opt_r.append(r)
            opt_ang.append(abs(ang))
            opt_err.append(np.rad2deg(max(le - vlo, vhi - he)))
    return {
        "n_geom": n_geom, "n_opt": len(opt_err), "errs": np.array(opt_err),
        "mean_lever_all": float(np.mean(all_r)),
        "mean_lever_opt": float(np.mean(opt_r)) if opt_r else float("nan"),
        "mean_ang_all": float(np.mean(all_ang)),
        "mean_ang_opt": float(np.mean(opt_ang)) if opt_ang else float("nan"),
    }


def report(label, res):
    n_opt, n_geom = res["n_opt"], res["n_geom"]
    print(f"{label}: {n_opt}/{n_geom} ({100*n_opt/n_geom:.2f}%) optimistic")
    if n_opt:
        print(f"  median error {np.median(res['errs']):.3f} deg, "
              f"max {res['errs'].max():.3f} deg")
        print(f"  mean lever arm: optimistic {res['mean_lever_opt']:.3f} m "
              f"vs all-sampled {res['mean_lever_all']:.3f} m")
        print(f"  mean |off-axis angle|: optimistic {res['mean_ang_opt']:.3f} rad "
              f"vs all-sampled {res['mean_ang_all']:.3f} rad")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mult", type=float, default=3.0,
                     help="delta half-width multiplier for the wider-box "
                          "check, relative to the paper scenario's +/-15 mm "
                          "(default 3, matching the paper text)")
    ap.add_argument("--n", type=int, default=2000)
    args = ap.parse_args()

    report(f"Baseline (delta = +/-{DELTA0*1000:.0f} mm)",
           optimism_sweep(make_box(1.0), args.n))
    report(f"delta x{args.mult:.1f} (+/-{DELTA0*args.mult*1000:.0f} mm)",
           optimism_sweep(make_box(args.mult), args.n))


if __name__ == "__main__":
    main()
