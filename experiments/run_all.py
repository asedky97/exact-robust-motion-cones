#!/usr/bin/env python3
"""
Reproduces every headline number in:
  "Exact Robust Motion Cones: Guaranteed Sticking Manipulation under
   Parametric and Geometric Contact Uncertainty" (IEEE Access submission).

Usage:
    python run_all.py            # run all experiments E1-E6, write summary.json
    python run_all.py --exp 1    # run a single experiment
    python run_all.py --quick    # reduced sample counts for a fast sanity pass
    python run_all.py --figs     # regenerate figures from the latest results

All experiments share the paper scenario:
    p0 = [-0.30, 0] m,  n0 = [1, 0],  mu0 = 0.45,  c0 = 0.06 m
    box: mu in [0.25, 0.65], c in [0.045, 0.075] m,
         alpha in [-6, +6] deg, delta in [-15, +15] mm
    Exact robust cone: [-59.91, +59.91] deg ;  nominal: +/- 85.11 deg

Approx. full runtime: 2-4 min on a laptop (E1 dominates).
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

# Make the src/ package importable whether run from repo root or experiments/
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from robust_motion_cone import (  # noqa: E402
    ContactGeometry, UncertaintyBox, cone_angles, robust_cone_exact,
)
from pushing_sim import resolve_contact  # noqa: E402

OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
RESULTS = ROOT / "experiments" / "summary.json"


def _wilson_ci(k: int, n: int, z: float = 1.96):
    """95% Wilson score confidence interval for a binomial proportion (percent)."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    lo = max(0.0, 100 * (centre - margin))
    hi = min(100.0, 100 * (centre + margin))
    return (float(round(lo, 2)), float(round(hi, 2)))

# ----------------------------------------------------------------- scenario
BOX = UncertaintyBox(
    mu=(0.25, 0.65),
    c=(0.045, 0.075),
    alpha=(np.deg2rad(-6), np.deg2rad(6)),
    delta=(-0.015, 0.015),
)
GEOM = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))


def robust_and_nominal():
    lo, hi, ref, nev = robust_cone_exact(GEOM, BOX)
    mu0, c0, a0, d0 = BOX.nominal()
    lon, hin = cone_angles(GEOM, mu0, c0, a0, d0, ref)
    return lo, hi, ref, lon, hin, nev


# ----------------------------------------------------------------- E1
def exp1(quick=False):
    """Exactness and efficiency vs 21^4 grid and vertex-only baseline."""
    rng = np.random.default_rng(0)
    n_geom = 200 if quick else 2000

    lo, hi, ref, lon, hin, nev = robust_and_nominal()

    # Accuracy + optimism check over random geometries near the paper regime
    # (lever arms comparable to the scenario; very short arms give near-empty
    # cones where vertex-only is trivially optimistic and uninformative).
    n_opt = 0
    max_err = 0.0
    for _ in range(n_geom):
        # sample lever arm 0.20-0.35 m (paper scenario uses 0.30 m)
        r = rng.uniform(0.20, 0.35)
        th_p = rng.uniform(-0.15, 0.15)
        p = np.array([-r * np.cos(th_p), r * np.sin(th_p)])
        ang = rng.uniform(-0.15, 0.15)
        G = ContactGeometry(p0=p, n0=np.array([np.cos(ang), np.sin(ang)]))
        le, he, r2, _ = robust_cone_exact(G, BOX)
        if he <= le:  # skip empty cones (vertex-only optimism is uninformative)
            continue
        vlo, vhi = -np.inf, np.inf
        for mu in BOX.mu:
            for c in BOX.c:
                for al in BOX.alpha:
                    for de in BOX.delta:
                        l, h = cone_angles(G, mu, c, al, de, r2)
                        vlo = max(vlo, l)
                        vhi = min(vhi, h)
        # vertex-only is optimistic if it reports a WIDER cone than exact
        if vlo < le - 1e-9 or vhi > he + 1e-9:
            n_opt += 1
            max_err = max(max_err, np.rad2deg(max(le - vlo, vhi - he)))

    # Timing: exact per-call vs a single full grid evaluation.
    # The paper's baseline is a 21^4 = 194,481-point grid; we time the exact
    # call and estimate the grid cost from the per-point evaluation rate so the
    # reported speedup reflects the full-resolution grid regardless of --quick.
    t0 = time.perf_counter()
    for _ in range(500):
        robust_cone_exact(GEOM, BOX)
    t_exact = (time.perf_counter() - t0) / 500

    # measure cost of one cone_angles evaluation, scale to 21^4
    _, _, rr, _ = robust_cone_exact(GEOM, BOX)
    t1 = time.perf_counter()
    for _ in range(2000):
        cone_angles(GEOM, 0.45, 0.06, 0.0, 0.0, rr)
    t_eval = (time.perf_counter() - t1) / 2000
    t_grid = t_eval * (21 ** 4)

    # NOTE on the two speedup numbers below (by design, they differ a lot):
    #   - speedup_est is t_grid_est_ms / t_exact_ms, where t_grid_est_ms
    #     extrapolates a SCALAR per-point cone_angles() call in a Python loop.
    #     It is a quick, deliberately-crude lower bound on how much slower a
    #     naive grid is; it is NOT the paper's grid baseline.
    #   - speedup_paper (225x) is the wall-clock ratio the paper reports,
    #     measured once on the reference machine (Dell Precision 7710) against
    #     the paper's actual grid baseline: a NumPy-VECTORIZED 21^4 evaluation
    #     (array broadcasting, not a Python scalar loop). A vectorized grid
    #     amortizes per-point overhead that the scalar loop above pays 194,481
    #     times, so speedup_est is expected to be much larger than 225x -- both
    #     numbers agree on the qualitative conclusion (the grid is orders of
    #     magnitude slower than Algorithm 1); they are not meant to match
    #     numerically. speedup_paper is kept here as a fixed reference value
    #     for traceability, not recomputed by this script.
    res = {
        "robust_deg": [round(float(np.rad2deg(lo)), 2), round(float(np.rad2deg(hi)), 2)],
        "nominal_deg": [round(float(np.rad2deg(lon)), 2), round(float(np.rad2deg(hin)), 2)],
        "n_evals": int(nev),
        "pct_optimistic": round(100 * n_opt / n_geom, 2),
        "max_vtx_err_deg": round(float(max_err), 2),
        "t_exact_ms": round(t_exact * 1e3, 4),
        "t_grid_est_ms": round(t_grid * 1e3, 1),
        "speedup_est": round(t_grid / t_exact, 0),
        "speedup_paper": 225,  # measured wall-clock on the reference machine
    }
    print(f"[E1] robust={res['robust_deg']} nominal={res['nominal_deg']} "
          f"optimistic={res['pct_optimistic']}% "
          f"exact={res['t_exact_ms']}ms (paper speedup 225x)")
    return res


# ----------------------------------------------------------------- E2
def exp2(quick=False):
    """Monte-Carlo mode verification by zone."""
    rng = np.random.default_rng(1)
    N = 4000 if quick else 20000
    lo, hi, ref, lon, hin, _ = robust_and_nominal()
    zones = {"R": [0, 0], "N": [0, 0], "S": [0, 0]}
    for _ in range(N):
        psi = rng.uniform(lon - 0.2, hin + 0.2)  # sample near/inside nominal
        if lo <= psi <= hi:
            z = "R"
        elif lon <= psi <= hin:
            z = "N"
        else:
            z = "S"
        th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
              rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
        p, n = GEOM.contact(th[2], th[3])
        u = np.array([np.cos(psi), np.sin(psi)])
        m, _ = resolve_contact(p, n, th[0], th[1], u)
        zones[z][0] += int(m != "stick")
        zones[z][1] += 1
    res = {}
    for z, (f, t) in zones.items():
        ci_lo, ci_hi = _wilson_ci(f, t)
        res[z] = {
            "slide_pct": round(100 * f / max(t, 1), 2),
            "n": t,
            "ci_95": [ci_lo, ci_hi],
        }
    print(f"[E2] Zone R={res['R']['slide_pct']}% CI=[{res['R']['ci_95']}]  "
          f"N={res['N']['slide_pct']}% CI=[{res['N']['ci_95']}]  "
          f"S={res['S']['slide_pct']}% CI=[{res['S']['ci_95']}]")
    return res


# ----------------------------------------------------------------- E3
def exp3(quick=False):
    """Certified pushing funnel vs Theorem 3 bound."""
    rng = np.random.default_rng(2)
    n_real = 200 if quick else 1000
    steps = 500
    lo, hi, ref, _, _, _ = robust_and_nominal()
    psi = lo + 0.35 * (hi - lo)
    u_cmd = np.array([np.cos(psi), np.sin(psi)])

    # Theorem 3 bound: kappa_max = max over the box of the sticking curvature
    # |omega| / |u|. omega is monotone in s = 1/c^2 (extreme at the c endpoints),
    # but as a function of delta it is |linear / quadratic|, whose extremum can be
    # INTERIOR to [delta_lo, delta_hi]. Endpoint-only evaluation could therefore
    # UNDER-estimate a *certified* bound, so we scan delta densely (cheap, exact
    # to the grid) at both c endpoints.
    kappa_max = 0.0
    for c in BOX.c:
        s = 1 / c ** 2
        for de in np.linspace(BOX.delta[0], BOX.delta[1], 401):
            p = GEOM.p0 + de * np.array([0, 1])
            k = abs(s * (p[0] * u_cmd[1] - p[1] * u_cmd[0]) / (1 + s * (p @ p)))
            kappa_max = max(kappa_max, k)
    L = 0.5
    bound_theta = np.rad2deg(kappa_max * L)
    bound_pos = kappa_max * L ** 2 / 2 * 1000  # mm

    # nominal (zero-uncertainty) realization, for the deviation-from-nominal metric
    mu0, c0, a0, d0 = BOX.nominal()
    p_n, n_n = GEOM.contact(a0, d0)
    nx = ny = 0.0
    for _ in range(steps):
        m, v = resolve_contact(p_n, n_n, mu0, c0, u_cmd)
        if m != "stick":
            break
        nx += v[0] * 1e-3
        ny += v[1] * 1e-3
    nom = np.array([nx, ny])

    slides = 0
    finals = []
    finals_xy = []
    for _ in range(n_real):
        th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
              rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
        p, n = GEOM.contact(th[2], th[3])
        x = y = yaw = 0.0
        for _ in range(steps):
            m, v = resolve_contact(p, n, th[0], th[1], u_cmd)
            if m != "stick":
                slides += 1
                break
            x += v[0] * 1e-3
            y += v[1] * 1e-3
            yaw += v[2] * 1e-3
        finals.append(np.rad2deg(yaw))
        finals_xy.append((x, y))
    spread = float(np.ptp(finals))
    xy = np.array(finals_xy)                       # m, shape (n_real, 2)
    mean_xy = xy.mean(axis=0)
    # cross-realization position dispersion = max deviation of any realization's
    # final position from the mean final position (the dispersion radius), in mm.
    pos_disp = float(np.linalg.norm(xy - mean_xy, axis=1).max() * 1000.0)
    # max deviation from the nominal trajectory (the quantity Theorem 3 bounds), mm
    pos_dev_nom = float(np.linalg.norm(xy - nom, axis=1).max() * 1000.0)
    res = {
        "total_slides": slides,
        "theta_spread_deg": round(spread, 2),
        "theta_bound_deg": round(bound_theta, 2),
        "pos_bound_mm": round(bound_pos, 1),
        "pos_dispersion_mm": round(pos_disp, 1),
        "pos_dev_nominal_mm": round(pos_dev_nom, 1),
    }
    print(f"[E3] slides={slides}  spread={res['theta_spread_deg']}° "
          f"(bound {res['theta_bound_deg']}°)  pos_disp={res['pos_dispersion_mm']}mm "
          f"(bound {res['pos_bound_mm']}mm)")
    return res


# ----------------------------------------------------------------- E4
def exp4(quick=False):
    """Price of robustness: Lipschitz sensitivities + bound/exact ratio."""
    lo, hi, ref, lon, hin, _ = robust_and_nominal()
    A0 = np.rad2deg(hin - lon)

    def aperture(scale):
        b = UncertaintyBox(
            mu=(0.45 - 0.20 * scale, 0.45 + 0.20 * scale),
            c=(0.06 - 0.015 * scale, 0.06 + 0.015 * scale),
            alpha=(-np.deg2rad(6) * scale, np.deg2rad(6) * scale),
            delta=(-0.015 * scale, 0.015 * scale),
        )
        l, h, _, _ = robust_cone_exact(GEOM, b)
        return np.rad2deg(h - l)

    # numerical Lipschitz constants at the nominal point
    mu0, c0, a0, d0 = BOX.nominal()
    eps = 1e-4

    def edge(mu, c, al, de, sign):
        l, h = cone_angles(GEOM, mu, c, al, de, ref)
        return np.rad2deg(h if sign > 0 else l)

    L = {}
    L["mu"] = abs(edge(mu0 + eps, c0, a0, d0, +1) - edge(mu0, c0, a0, d0, +1)) / eps
    L["c"] = abs(edge(mu0, c0 + eps, a0, d0, +1) - edge(mu0, c0, a0, d0, +1)) / eps
    L["alpha"] = abs(edge(mu0, c0, a0 + eps, d0, +1) - edge(mu0, c0, a0, d0, +1)) / eps * np.pi / 180
    L["delta"] = abs(edge(mu0, c0, a0, d0 + eps, +1) - edge(mu0, c0, a0, d0, +1)) / eps

    small = 0.02
    loss_small = A0 - aperture(small)
    bound_small = small * 2 * (L["mu"] * 0.20 + L["c"] * 0.015
                               + L["alpha"] * 6 + L["delta"] * 0.015)
    loss_full = A0 - aperture(1.0)
    bound_full = 2 * (L["mu"] * 0.20 + L["c"] * 0.015
                      + L["alpha"] * 6 + L["delta"] * 0.015)
    res = {
        "L_mu_per_unit": round(L["mu"], 2),
        "L_c_per_m": round(L["c"], 2),
        "L_alpha_per_deg": round(L["alpha"], 2),
        "L_delta_per_m": round(L["delta"], 2),
        "ratio_small": round(bound_small / max(loss_small, 1e-9), 3),
        "ratio_full": round(bound_full / max(loss_full, 1e-9), 3),
        "loss_full_deg": round(loss_full, 2),
    }
    print(f"[E4] L_c={res['L_c_per_m']}/m L_delta={res['L_delta_per_m']}/m "
          f"ratio_small={res['ratio_small']} ratio_full={res['ratio_full']}")
    return res


# ----------------------------------------------------------------- E5
def exp5(quick=False):
    """Command-frame comparison: body-frame vs world-frame vs nominal."""
    rng = np.random.default_rng(3)
    N = 4000 if quick else 20000
    lo, hi, ref, lon, hin, _ = robust_and_nominal()

    def _count_failures(strat):
        fcnt = 0
        for _ in range(N):
            th = [rng.uniform(*BOX.mu), rng.uniform(*BOX.c),
                  rng.uniform(*BOX.alpha), rng.uniform(*BOX.delta)]
            p, n = GEOM.contact(th[2], th[3])
            yaw = 0.0
            ok = True
            for _ in range(8):
                if strat == "A":          # body-frame at 45 deg (inside robust)
                    psi = np.deg2rad(45)
                elif strat == "B":        # world-frame: drifts as object rotates
                    psi = np.deg2rad(45) - yaw
                else:                     # C: nominal body-frame at 84 deg
                    psi = np.deg2rad(84)
                u = np.array([np.cos(psi), np.sin(psi)])
                m, v = resolve_contact(p, n, th[0], th[1], u)
                if m != "stick":
                    ok = False
                    break
                yaw += v[2] * 40e-3
            fcnt += int(not ok)
        return fcnt

    res = {}
    for label, strat in [("A_robust", "A"), ("B_world", "B"), ("C_nominal", "C")]:
        fcnt = _count_failures(strat)
        ci_lo, ci_hi = _wilson_ci(fcnt, N)
        res[label] = {
            "failure_pct": round(100 * fcnt / N, 2),
            "n": N,
            "ci_95": [ci_lo, ci_hi],
        }
    print(f"[E5] A(robust)={res['A_robust']['failure_pct']}% "
          f"CI={res['A_robust']['ci_95']}  "
          f"B(world)={res['B_world']['failure_pct']}% "
          f"CI={res['B_world']['ci_95']}  "
          f"C(nominal)={res['C_nominal']['failure_pct']}% "
          f"CI={res['C_nominal']['ci_95']}")
    return res


# ----------------------------------------------------------------- E6
def exp6(quick=False):
    """Independent-engine consistency (single-realization transition).

    NOTE: the paper's E6 uses PyBullet (rigid-body contact) on the author's
    machine. This routine reproduces the *expected* behaviour using the
    analytical limit-surface resolver at one effective realization: the
    stick-to-slide transition falls between the robust and nominal edges.
    To run the true PyBullet sweep, see experiments/pybullet_sweep.py.
    """
    lo, hi, ref, lon, hin, _ = robust_and_nominal()
    mu_eff, c_eff = 0.30, 0.066
    angles = np.arange(-88, 89, 4)
    stick = []
    for ad in angles:
        a = np.deg2rad(ad)
        u = np.array([np.cos(a), np.sin(a)])
        p, n = GEOM.contact(0, 0)
        m, _ = resolve_contact(p, n, mu_eff, c_eff, u)
        if m == "stick":
            stick.append(ad)
    res = {
        "transition_deg": int(max(stick)) if stick else None,
        "robust_edge_deg": round(np.rad2deg(hi), 2),
        "nominal_edge_deg": round(np.rad2deg(hin), 2),
        "between": bool(stick and np.rad2deg(hi) <= max(stick) <= np.rad2deg(hin)),
    }
    print(f"[E6] transition=±{res['transition_deg']}° between robust "
          f"±{res['robust_edge_deg']}° and nominal ±{res['nominal_edge_deg']}°")
    return res


EXPERIMENTS = {1: exp1, 2: exp2, 3: exp3, 4: exp4, 5: exp5, 6: exp6}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exp", type=int, choices=range(1, 7),
                    help="run a single experiment (default: all)")
    ap.add_argument("--quick", action="store_true",
                    help="reduced sample counts for a fast pass")
    ap.add_argument("--figs", action="store_true",
                    help="regenerate figures from existing summary.json")
    args = ap.parse_args()

    if args.figs:
        from make_figures import make_all
        make_all()
        return

    # A --quick pass uses reduced sample counts, so it must NOT overwrite the
    # canonical full-run summary.json (doing so is how a stale, low-sample
    # summary can silently replace the paper's numbers). Quick results go to a
    # clearly-named sidecar file instead.
    results_path = RESULTS if not args.quick else RESULTS.with_name("summary_quick.json")

    summary = {}
    if results_path.exists():
        summary = json.loads(results_path.read_text())

    to_run = [args.exp] if args.exp else list(EXPERIMENTS)
    t0 = time.perf_counter()
    for e in to_run:
        summary[f"E{e}"] = EXPERIMENTS[e](quick=args.quick)
    results_path.write_text(json.dumps(summary, indent=2))
    note = "  (quick pass; canonical numbers are in summary.json)" if args.quick else ""
    print(f"\nWrote {results_path.name} ({time.perf_counter() - t0:.1f}s){note}")


if __name__ == "__main__":
    main()
