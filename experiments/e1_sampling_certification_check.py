#!/usr/bin/env python3
"""
Reproduces Sec. V.A's "sampling cannot certify worst-case safety" claim.

Phi_minus = max_Theta phi_minus(Theta) is a maximum over the uncertainty
box; Phi_plus = min_Theta phi_plus(Theta) is a minimum. A sample-based
estimate -- draw n realizations of Theta, take the sample max of
phi_minus and the sample min of phi_plus -- is therefore always at least
as wide as the true robust cone: a finite sample can only fail to reach
the true worst-case realization, never exceed it. The resulting interval
looks like a valid safe range from the sample alone, but is not
certified: a not-yet-sampled Theta could realize a stricter edge than any
edge seen so far, and closing that gap only gets slower as n grows
(bounded by extreme-value convergence, not the usual 1/n Monte-Carlo
rate for means).

Usage:
    python e1_sampling_certification_check.py [--seed 0]
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import ContactGeometry, UncertaintyBox, cone_angles, robust_cone_exact  # noqa

GEOM = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))
BOX = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                      alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))


def sample_estimate(n, rng, ref):
    phi_lo = np.empty(n)
    phi_hi = np.empty(n)
    for i in range(n):
        mu = rng.uniform(*BOX.mu)
        c = rng.uniform(*BOX.c)
        al = rng.uniform(*BOX.alpha)
        de = rng.uniform(*BOX.delta)
        l, h = cone_angles(GEOM, mu, c, al, de, ref)
        phi_lo[i] = l
        phi_hi[i] = h
    # sample-based estimate of Phi_minus (a max) and Phi_plus (a min)
    sample_lo = np.rad2deg(phi_lo.max())
    sample_hi = np.rad2deg(phi_hi.min())
    return sample_lo, sample_hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    true_lo, true_hi, ref, _ = robust_cone_exact(GEOM, BOX)
    true_lo, true_hi = np.rad2deg(true_lo), np.rad2deg(true_hi)
    true_aperture = true_hi - true_lo
    print(f"True robust cone: [{true_lo:.2f}, {true_hi:.2f}] deg "
          f"(aperture {true_aperture:.2f} deg)")

    rng = np.random.default_rng(args.seed)
    for n in (5_000, 20_000):
        sample_lo, sample_hi = sample_estimate(n, rng, ref)
        excess = (sample_hi - sample_lo) - true_aperture
        print(f"n={n:>6}: sample=[{sample_lo:.2f}, {sample_hi:.2f}] deg  "
              f"excess aperture over true = {excess:.2f} deg")


if __name__ == "__main__":
    main()
