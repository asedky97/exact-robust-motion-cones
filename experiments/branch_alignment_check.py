#!/usr/bin/env python3
"""
Branch-alignment and quadratic-degeneracy diagnostics over 10,000 random
geometries -- reproduces two specific supporting empirical claims that are
NOT covered by run_all.py's six headline experiments (E1-E6):

  1. Remark 1 (Section II.C): how often does the initial reference suffice
     for branch alignment (Corollary 1), and how often is one re-centre
     needed? Result: r0 alone has sufficed in every sampled geometry so far
     (both the realistic default sampling and --wide adversarial sampling) --
     the recenter fallback itself has never been triggered and so remains
     empirically untested, not just unproven (see proofs_supplement.md).
  2. The numerical-implementation note in Section III.C (Theorem 2): across
     how many geometries does the delta-critical-point quadratic degenerate
     (near-zero leading coefficient)?

Recentre count is inferred from n_evals: each internal _eval_candidates
pass does at most 16 evaluations (Proposition 1), and robust_cone_exact
accumulates n_evals additively across up to 3 passes, so
  n_evals <= 16   -> 1 pass  (0 recentres)
  16 < n_evals <= 32 -> 2 passes (1 recentre)
  n_evals > 32    -> 3 passes (2 recentres)
This reads a value robust_cone_exact already returns; it does not modify or
duplicate any tested logic in robust_motion_cone.py.

The degeneracy check re-derives the same quadratic leading coefficient
Theorem 2 defines (Bε x Cε), using the same public helper functions
(perp, rot, friction_edges, cross2) that _delta_critical_points itself is
built from, for each of the up to 4 (edge, c) candidate quadratics per
geometry -- matching Algorithm 1's own candidate enumeration.

Usage:
    python branch_alignment_check.py [--n 10000]
    python branch_alignment_check.py --wide [--n 10000]   # adversarial sweep
                                                            # (lever arms down to
                                                            # 0.03 m, normal
                                                            # angles spanning the
                                                            # full circle) cited
                                                            # in Remark 1
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import (  # noqa
    ContactGeometry, UncertaintyBox, robust_cone_exact,
    perp, rot, friction_edges, cross2,
)

BOX = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                     alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))


def quadratic_leading_coeff(geom, mu_lo, c, alpha_val, which):
    """B x C for the edge quadratic (Theorem 2, Eq. 7-8), computed exactly as
    _delta_critical_points does internally (same helper functions)."""
    s = 1.0 / c ** 2
    n = rot(alpha_val) @ geom.n0
    fm, fp = friction_edges(n, mu_lo)
    f = fm if which == "minus" else fp
    q1, q2 = perp(geom.p0), perp(geom.t0)
    k0, k1 = float(q1 @ f), float(q2 @ f)
    B = s * (k1 * q1 + k0 * q2)
    C = s * k1 * q2
    return cross2(B, C)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10_000)
    ap.add_argument("--wide", action="store_true",
                     help="adversarial sampling (lever arms down to 0.03 m, "
                          "normal angles spanning the full circle) instead of "
                          "the paper-scenario regime")
    args = ap.parse_args()
    N = args.n

    rng = np.random.default_rng(99 if args.wide else 4)
    pass_counts = {1: 0, 2: 0, 3: 0}
    geoms_with_degenerate = 0
    total_degenerate_quadratics = 0
    total_quadratics = 0

    for _ in range(N):
        if args.wide:
            r = rng.uniform(0.03, 0.5)
            th_p = rng.uniform(-3.0, 3.0)
            ang = rng.uniform(-3.0, 3.0)
        else:
            r = rng.uniform(0.20, 0.35)
            th_p = rng.uniform(-0.15, 0.15)
            ang = rng.uniform(-0.15, 0.15)
        p = np.array([-r * np.cos(th_p), r * np.sin(th_p)])
        geom = ContactGeometry(p0=p, n0=np.array([np.cos(ang), np.sin(ang)]))

        _, _, _, nev = robust_cone_exact(geom, BOX)
        n_pass = 1 if nev <= 16 else (2 if nev <= 32 else 3)
        pass_counts[n_pass] += 1

        geom_degenerate = False
        for which, alpha_val in (("minus", BOX.alpha[1]), ("plus", BOX.alpha[0])):
            for c in set(BOX.c):
                total_quadratics += 1
                leading = quadratic_leading_coeff(geom, BOX.mu[0], c, alpha_val, which)
                if abs(leading) < 1e-14:
                    total_degenerate_quadratics += 1
                    geom_degenerate = True
        if geom_degenerate:
            geoms_with_degenerate += 1

    print(f"Geometries: {N} ({'adversarial/wide' if args.wide else 'paper-scenario'} sampling)")
    print(f"Recentre passes: 1st-try {pass_counts[1]} ({100*pass_counts[1]/N:.2f}%), "
          f"1 recentre {pass_counts[2]} ({100*pass_counts[2]/N:.2f}%), "
          f"2 recentres {pass_counts[3]} ({100*pass_counts[3]/N:.2f}%)")
    print(f"Max recentres observed: {max(k for k, v in pass_counts.items() if v > 0) - 1}")
    print(f"Degenerate quadratics: {total_degenerate_quadratics} / {total_quadratics} "
          f"candidate evaluations, across {geoms_with_degenerate} / {N} geometries")


if __name__ == "__main__":
    main()
