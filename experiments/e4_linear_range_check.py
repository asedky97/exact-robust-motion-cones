#!/usr/bin/env python3
"""
Reproduces the "for lambda <~ 0.1, the linear estimate approximates the true
loss to within ~5%" claim (Remark 4 and Section V.D / E4's text) across a
range of box scales -- run_all.py::exp4 only evaluates two points, lambda
= 0.02 ("small") and lambda = 1.0 ("full"), not the claimed lambda <= 0.1
regime.

Replicates exp4's exact aperture() and Lipschitz-bound computation and
sweeps the box scale from max-lambda/steps up to max-lambda.

Usage:
    python e4_linear_range_check.py [--max-lambda 0.1] [--steps 5]
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

GEOM = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))
BOX = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                     alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))


def aperture(scale):
    """Exact aperture at box scale `scale` (identical to run_all.py::exp4)."""
    b = UncertaintyBox(
        mu=(0.45 - 0.20 * scale, 0.45 + 0.20 * scale),
        c=(0.06 - 0.015 * scale, 0.06 + 0.015 * scale),
        alpha=(-np.deg2rad(6) * scale, np.deg2rad(6) * scale),
        delta=(-0.015 * scale, 0.015 * scale),
    )
    l, h, _, _ = robust_cone_exact(GEOM, b)
    return np.rad2deg(h - l)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-lambda", type=float, default=0.10,
                     help="largest box scale to test (default 0.10, matching "
                          "Remark 4's claimed regime)")
    ap.add_argument("--steps", type=int, default=5)
    args = ap.parse_args()

    lo, hi, ref, nev = robust_cone_exact(GEOM, BOX)
    mu0, c0, a0, d0 = BOX.nominal()
    lon, hin = cone_angles(GEOM, mu0, c0, a0, d0, ref)
    A0 = np.rad2deg(hin - lon)
    eps = 1e-4

    def edge(mu, c, al, de, sign):
        l, h = cone_angles(GEOM, mu, c, al, de, ref)
        return np.rad2deg(h if sign > 0 else l)

    # Identical Lipschitz-constant computation to run_all.py::exp4.
    L = {}
    L["mu"] = abs(edge(mu0 + eps, c0, a0, d0, +1) - edge(mu0, c0, a0, d0, +1)) / eps
    L["c"] = abs(edge(mu0, c0 + eps, a0, d0, +1) - edge(mu0, c0, a0, d0, +1)) / eps
    L["alpha"] = (abs(edge(mu0, c0, a0 + eps, d0, +1) - edge(mu0, c0, a0, d0, +1))
                  / eps * np.pi / 180)
    L["delta"] = abs(edge(mu0, c0, a0, d0 + eps, +1) - edge(mu0, c0, a0, d0, +1)) / eps

    def bound(scale):
        return scale * 2 * (L["mu"] * 0.20 + L["c"] * 0.015
                            + L["alpha"] * 6 + L["delta"] * 0.015)

    print(f"{'lambda':>8} {'loss(deg)':>10} {'bound(deg)':>11} {'rel.err':>8}")
    worst = 0.0
    for scale in np.linspace(args.max_lambda / args.steps, args.max_lambda, args.steps):
        loss = A0 - aperture(scale)
        b = bound(scale)
        err_pct = 100 * abs(b - loss) / loss
        worst = max(worst, err_pct)
        print(f"{scale:8.3f} {loss:10.3f} {b:11.3f} {err_pct:7.2f}%")
    print(f"\nWorst-case relative error up to lambda={args.max_lambda}: {worst:.2f}%")


if __name__ == "__main__":
    main()
