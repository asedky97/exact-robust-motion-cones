#!/usr/bin/env python3
"""
Reproduces the "50%-failure horizon at ~44 deg of cumulative rotation" claim
for Strategy B (Section V.E) that run_all.py::exp5 does not expose -- exp5
only counts whether each trial fails anywhere across its 8 pushes, not the
cumulative body-frame rotation at the moment of failure.

Replicates exp5's exact methodology, including running strategies in the
same order (A, then B, then C) on the SAME numpy Generator instance, since
exp5 does not reset the RNG between strategies -- testing strategy B in
isolation with a fresh seed draws from the wrong part of the random stream
and silently gives a different (wrong) answer.

For each of strategy B's failing trials, records the cumulative |yaw| at
the step where it first fails to stick; the median of that distribution is
the reported "50%-failure horizon".

Usage:
    python e5_rotation_horizon_check.py [--n 20000]
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import ContactGeometry, UncertaintyBox  # noqa
from pushing_sim import resolve_contact  # noqa

GEOM = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))
BOX = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                     alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))


def run_strategy(rng, strat, n, track_horizon=False):
    """Identical mechanics to run_all.py::exp5's _count_failures, optionally
    recording the cumulative |yaw| (deg) at the failing step of each trial
    that fails."""
    fcnt = 0
    fail_yaws = []
    for _ in range(n):
        th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
              rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
        p, n_ = GEOM.contact(th[2], th[3])
        yaw = 0.0
        ok = True
        for _ in range(8):
            if strat == "A":
                psi = np.deg2rad(45)
            elif strat == "B":
                psi = np.deg2rad(45) - yaw
            else:
                psi = np.deg2rad(84)
            u = np.array([np.cos(psi), np.sin(psi)])
            m, v = resolve_contact(p, n_, th[0], th[1], u)
            if m != "stick":
                ok = False
                if track_horizon:
                    fail_yaws.append(np.rad2deg(abs(yaw)))
                break
            yaw += v[2] * 40e-3
        fcnt += int(not ok)
    return fcnt, fail_yaws


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20_000)
    args = ap.parse_args()
    N = args.n

    rng = np.random.default_rng(3)  # same seed and same Generator instance
    fA, _ = run_strategy(rng, "A", N)  # must run first: consumes rng state
    print(f"A fails: {fA}/{N} ({100*fA/N:.2f}%)")

    fB, fail_yaws = run_strategy(rng, "B", N, track_horizon=True)
    fy = np.array(fail_yaws)
    print(f"B fails: {fB}/{N} ({100*fB/N:.2f}%)")
    print(f"  50%-failure horizon (median cumulative |yaw| at failure): "
          f"{np.median(fy):.2f} deg")
    print(f"  25th/75th percentile: {np.percentile(fy, 25):.2f} / "
          f"{np.percentile(fy, 75):.2f} deg")

    fC, _ = run_strategy(rng, "C", N)  # must run to match exp5's rng sequence
    print(f"C fails: {fC}/{N} ({100*fC/N:.2f}%)")


if __name__ == "__main__":
    main()
