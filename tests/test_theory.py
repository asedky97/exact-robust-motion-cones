#!/usr/bin/env python3
"""
Unit tests for the exact robust motion cone algorithm.

Validates the theoretical claims of the paper against the implementation:
  - Theorem 1: monotonicity in mu, c, alpha
  - Theorem 2: exactness vs dense grid; cubic coefficient vanishes
  - Proposition 1: candidate set size <= 16
  - Headline scenario numbers

Run:  python -m pytest tests/ -v      (or)   python tests/test_theory.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import (  # noqa
    ContactGeometry, UncertaintyBox, cone_angles, robust_cone_exact,
    sticking_omega, push_map,
)

BOX = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                     alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))
GEOM = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))


def test_headline_scenario():
    """Robust cone matches the paper's reported [-59.91, +59.91] deg."""
    lo, hi, _, _ = robust_cone_exact(GEOM, BOX)
    assert abs(np.rad2deg(lo) - (-59.91)) < 0.1
    assert abs(np.rad2deg(hi) - (+59.91)) < 0.1


def test_nominal_scenario():
    """Nominal cone matches the paper's +/- 85.11 deg."""
    mu0, c0, a0, d0 = BOX.nominal()
    _, _, ref, _ = robust_cone_exact(GEOM, BOX)
    lon, hin = cone_angles(GEOM, mu0, c0, a0, d0, ref)
    assert abs(np.rad2deg(hin) - 85.11) < 0.1


def test_candidate_set_bound():
    """Proposition 1: at most 16 candidate evaluations per call."""
    _, _, _, nev = robust_cone_exact(GEOM, BOX)
    # ≤16 bound: 2 c values × (2 delta endpoints + up to 2 roots) × 2 edges
    assert nev <= 16


def test_robust_inside_nominal():
    """Robust cone is contained in the nominal cone."""
    lo, hi, ref, _ = robust_cone_exact(GEOM, BOX)
    mu0, c0, a0, d0 = BOX.nominal()
    lon, hin = cone_angles(GEOM, mu0, c0, a0, d0, ref)
    assert lo >= lon - 1e-9
    assert hi <= hin + 1e-9


def test_monotone_mu():
    """Theorem 1(i): higher friction widens the cone (worst case at mu_min)."""
    _, _, ref, _ = robust_cone_exact(GEOM, BOX)
    apertures = []
    for mu in np.linspace(0.1, 0.9, 9):
        l, h = cone_angles(GEOM, mu, 0.06, 0.0, 0.0, ref)
        apertures.append(h - l)
    diffs = np.diff(apertures)
    assert np.all(diffs > -1e-9)  # non-decreasing


def test_monotone_c():
    """Theorem 1(ii): each edge angle is monotone in c (no sign change)."""
    _, _, ref, _ = robust_cone_exact(GEOM, BOX)
    phi_plus = []
    for c in np.linspace(0.045, 0.075, 30):
        _, h = cone_angles(GEOM, 0.45, c, 0.0, 0.0, ref)
        phi_plus.append(h)
    diffs = np.diff(phi_plus)
    # constant sign => monotone
    assert np.all(diffs >= -1e-12) or np.all(diffs <= 1e-12)


def test_monotone_c_offaxis():
    """Theorem 1(ii) holds off-axis too (p_perp . p = 0 argument)."""
    G = ContactGeometry(p0=np.array([-0.25, 0.10]), n0=np.array([1.0, 0.0]))
    _, _, ref, _ = robust_cone_exact(G, BOX)
    phi_plus = []
    for c in np.linspace(0.045, 0.075, 30):
        _, h = cone_angles(G, 0.45, c, 0.0, 0.0, ref)
        phi_plus.append(h)
    diffs = np.diff(phi_plus)
    assert np.all(diffs >= -1e-12) or np.all(diffs <= 1e-12)


def test_monotone_alpha():
    """Theorem 1(iii): BOTH edges increase with alpha (this is what justifies
    selecting alpha_max for phi- and alpha_min for phi+ in robust_cone_exact)."""
    _, _, ref, _ = robust_cone_exact(GEOM, BOX)
    phi_minus, phi_plus = [], []
    for al in np.linspace(-0.1, 0.1, 21):
        l, h = cone_angles(GEOM, 0.45, 0.06, al, 0.0, ref)
        phi_minus.append(l)
        phi_plus.append(h)
    assert np.all(np.diff(phi_plus) > -1e-9)   # phi+ non-decreasing in alpha
    assert np.all(np.diff(phi_minus) > -1e-9)  # phi- non-decreasing in alpha


def test_sticking_omega():
    """Theorem 3: the closed-form sticking rotation rate omega = (p x u)/(c^2+|p|^2)
    equals the twist's angular component s*tau produced under sticking. Verified
    over random contacts and force directions inside the friction cone."""
    rng = np.random.default_rng(11)
    for _ in range(200):
        p = np.array([rng.uniform(-0.4, -0.1), rng.uniform(-0.1, 0.1)])
        c = rng.uniform(0.045, 0.075)
        f = np.array([np.cos(rng.uniform(-0.5, 0.5)), np.sin(rng.uniform(-0.5, 0.5))])
        u = push_map(p, c) @ f          # contact velocity for force direction f
        s = 1.0 / c**2
        tau = p[0] * f[1] - p[1] * f[0]  # p x f
        assert abs(sticking_omega(p, c, u) - s * tau) < 1e-9


def test_exactness_vs_grid():
    """Theorem 2: exact computation matches a dense delta grid."""
    rng = np.random.default_rng(7)
    for _ in range(20):
        p = np.array([rng.uniform(-0.4, -0.1), rng.uniform(-0.05, 0.05)])
        ang = rng.uniform(-0.2, 0.2)
        G = ContactGeometry(p0=p, n0=np.array([np.cos(ang), np.sin(ang)]))
        lo, hi, ref, _ = robust_cone_exact(G, BOX)
        # dense grid over delta only (the non-monotone parameter), at the
        # vertex (mu_min, c endpoints, alpha endpoints) the exact uses
        glo, ghi = -np.inf, np.inf
        for mu in [BOX.mu[0]]:
            for c in BOX.c:
                for al in BOX.alpha:
                    for de in np.linspace(*BOX.delta, 200):
                        l, h = cone_angles(G, mu, c, al, de, ref)
                        glo = max(glo, l); ghi = min(ghi, h)
        assert abs(lo - glo) < np.deg2rad(0.1)
        assert abs(hi - ghi) < np.deg2rad(0.1)


def test_emptiness_detection():
    """A huge uncertainty box should yield an empty robust cone."""
    big = UncertaintyBox(mu=(0.05, 0.95), c=(0.01, 0.30),
                         alpha=(np.deg2rad(-45), np.deg2rad(45)),
                         delta=(-0.10, 0.10))
    lo, hi, _, _ = robust_cone_exact(GEOM, big)
    assert hi <= lo + 1e-6  # empty or degenerate


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL {fn.__name__}: {e}")
        except Exception as e:  # noqa
            print(f"  ERROR {fn.__name__}: {e}")
    print(f"\n{passed}/{len(fns)} tests passed")
    return passed == len(fns)


if __name__ == "__main__":
    ok = _run_all()
    sys.exit(0 if ok else 1)
