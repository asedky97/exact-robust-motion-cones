# Reproducibility

This document gives the exact steps to reproduce every numerical result and
figure in *Exact Robust Motion Cones* (IEEE Access submission).

## Environment

| Component | Version used | Minimum |
|-----------|--------------|---------|
| Python    | 3.11         | 3.9     |
| NumPy     | 1.26         | 1.21    |
| matplotlib| 3.8          | 3.5     |
| PyBullet (E6 only) | 3.2.7 | 3.2 |
| OS        | Windows 11 / Ubuntu 22.04 | any |

No GPU, no special hardware, and no internet access are required. The core
algorithm depends only on NumPy. Total wall-clock time for the full campaign is
2-4 minutes on a 2020-era laptop (Experiment 1 dominates).

## One-command reproduction

```bash
pip install -r requirements.txt
cd experiments
python run_all.py          # writes summary.json with every headline number
python make_figures.py     # writes all figures to ../figures/
```

To reproduce the independent-engine experiment (E6):

```bash
pip install pybullet
python experiments/pybullet_sweep.py
```

## Determinism

- Every Monte-Carlo routine seeds `numpy.random.default_rng(k)` with a fixed
  integer `k` (shown in `run_all.py`), so results are bit-reproducible on a
  given NumPy version.
- Across NumPy versions, sampling-based numbers (E2, E3, E5) may vary by
  < 0.5 percentage points; the exact-algorithm numbers (E1 cone edges, E4
  Lipschitz constants) are deterministic to machine precision.
- The headline cone `[-59.91, +59.91]` deg and nominal `+/- 85.11` deg are
  computed in closed form and are exact.
- **Not deterministic, by design:** `summary.json`'s wall-clock timing fields
  (`E1.t_exact_ms`, `E1.t_grid_est_ms`, `E1.speedup_est`) are a live
  microbenchmark and will differ between machines and even between runs on
  the same machine; this is expected, not a reproduction failure.
  `speedup_paper` (225x) is the one timing-derived number that is fixed text
  in the paper itself, not something `run_all.py` re-measures; see the
  `speedup_est` vs. `speedup_paper` note further below.

## Mapping results to the paper

| Paper item | Produced by | Key in summary.json |
|------------|-------------|---------------------|
| Sec. IV worked example | `run_all.py --exp 1` | `E1.robust_deg`, `E1.nominal_deg` |
| Table (speedup) | `run_all.py --exp 1` | `E1.speedup_est`, `E1.t_exact_ms` |
| Fig. (sensitivity) | `make_figures.py` | N/A (`fig_sensitivity`) |
| E2 zone failure | `run_all.py --exp 2` | `E2.R/N/S` |
| E3 funnel | `run_all.py --exp 3` | `E3.total_slides`, `E3.theta_spread_deg` |
| E4 Lipschitz + ratio | `run_all.py --exp 4` | `E4.L_*`, `E4.ratio_small/full` |
| E5 command-frame | `run_all.py --exp 5` | `E5.A/B/C` |
| E6 engine transition | `pybullet_sweep.py` | (PyBullet) / `E6` (analytical) |
| E6(iii) box-spanning independent-solver check | `pybullet_box_sweep.py` | N/A (`pybullet_box_sweep_results.json`) |
| Fig. 7 (independent-engine composite) | see "Reproducing Fig. 7" below | N/A |
| Remark 1 (branch alignment) / §III.C degeneracy note | `branch_alignment_check.py` (add `--wide` for the adversarial sweep also cited in Remark 1) | N/A (console output) |
| E1 median error, lever-arm/off-axis correlation, wider-δ-box sensitivity | `e1_optimism_diagnostics.py` | N/A (console output) |
| §V.A sampling-cannot-certify claim (sample max/min vs. true robust cone) | `e1_sampling_certification_check.py` | N/A (console output) |
| E5 rotation-failure horizon (~44°) | `e5_rotation_horizon_check.py` | N/A (console output) |
| Remark 4 / E4 linear-estimate accuracy across λ ≤ 0.1 | `e4_linear_range_check.py` | N/A (console output) |

## Reproducing Fig. 7 (independent-engine validation composite)

Fig. 7 combines two Simscape validations (stick/slide sweep, mu-sweep) into
one 2-panel figure. Each underlying sweep script exports its results as
JSON; `make_fig7_composite.py` reads both and renders the composite in the
paper's matplotlib style:

```bash
# In MATLAB (Simscape + Simulink required), from supplementary/:
robust_motion_cone_simscape   # writes supplementary/main_sweep_results.json
cone_simscape_sweep           # writes supplementary/mu_sweep_results.json

# Then, in Python, from code/experiments/:
python make_fig7_composite.py             # reads both JSON files by default,
                                           # writes figures/fig5_engine_sweep_validated.png/.pdf
```

Expected: 84.1% (37/44) Simscape stick/slide accuracy, 88.9% precision, 85.7%
recall; Simscape mu-sweep max error 0.00 deg.

**A PyBullet panel is not part of Fig. 7.** An identification approach that
fits effective per-realization parameters (mu_eff, c_eff) and reports a
predicted-vs-observed residual is described, along with why the residual is
not used as validation evidence, in
`supplementary/pybullet_identification_retraction.md`. `python
experiments/pybullet_sweep.py` still runs and still writes
`pybullet_results.json`; its fitted-parameter output is not treated as
validated and is not part of Fig. 7.

## Reproducing E6(iii): box-spanning independent-solver check

Unlike `pybullet_sweep.py` (one fitted realization), `pybullet_box_sweep.py`
tests the robust cone itself: mu is sampled uniformly across the full box
`[0.25, 0.65]` via `changeDynamics`, with commands drawn from the same three
zones as E2 (R = inside robust, N = inside nominal but outside robust,
S = outside nominal).

```bash
python experiments/pybullet_box_sweep.py          # 300 trials/zone (paper regime)
python experiments/pybullet_box_sweep.py --n 100   # faster spot-check
```

Expected: Zone R 0.00% slide (95% CI [0.00, 1.26]%), Zone N 29.33% slide,
Zone S 100.00% slide. This test uses a different contact geometry from the
"paper scenario" used everywhere else (robust cone [-6.70, +6.70] deg,
nominal +/-60.62 deg for this geometry, not +/-59.91/+/-85.11); see the
board-sizing note below for why.

Only mu is swept: `alpha` and `delta` are held at nominal, and `c` is not
independently controllable in PyBullet's rigid-body primitive, so it is
fixed at one value rather than varied. That value is not arbitrary:
`pybullet_box_sweep_world.py`'s `BOARD_HALF` is chosen so the board's
footprint gives an implied `c` (via the same radius-of-gyration formula the
paper uses for a uniform rectangular support, Sec. II.A) that falls inside
the assumed c-box, and the contact point (`CONTACT_X`) sits at that board's
edge to match. This matters: an earlier version of this script used a much
larger board whose implied `c` was roughly 2-3x the assumed box, which
produced spurious slides inside Zone R that looked like a refutation of the
robust-cone guarantee but were actually a mismatch between the test
geometry and the model it was supposed to be testing. Fitting the observed
slide/stick pattern against a swept range of `c` values confirmed this:
agreement between the analytical prediction and PyBullet peaked at 98.3%
around `c ≈ 0.25` m (the actual implied `c` of that earlier board), far
outside the assumed box, and fell off on either side of that peak: a clean
signature of a geometry mismatch, not a modeling error. `BOARD_HALF` now
keeps the same aspect ratio, scaled down so its implied `c` sits inside the
box the test is meant to check.

## Reproducing Remark 1 (branch alignment) and the quadratic-degeneracy note

```bash
python experiments/branch_alignment_check.py           # 10,000 geometries, paper regime
python experiments/branch_alignment_check.py --wide     # 10,000 geometries, adversarial regime (full-circle normal angle, lever arms down to 3 cm)
python experiments/branch_alignment_check.py --n 1000   # faster spot-check
```

Reports how many of the sampled geometries required a re-centre (Corollary 1)
and how many candidate quadratics (Theorem 2) were degenerate. Expected: 0
recentres and 0 degenerate quadratics in both the paper regime and the
adversarial regime, meaning the one-recenter fallback itself is not
exercised by either sweep; see `supplementary/proofs_supplement.md` for what
that does and doesn't establish about Corollary 1.

## Reproducing E1's baseline-comparison sub-claims

`run_all.py::exp1` only reports the two headline E1 numbers (1.85% optimistic,
0.24 deg max error). The paragraph's other sub-claims -- the median error, the
correlation of optimism with lever arm / off-axis angle, and the wider-delta-box
sensitivity -- are reproduced separately:

```bash
python experiments/e1_optimism_diagnostics.py             # delta x3 (paper text)
python experiments/e1_optimism_diagnostics.py --mult 2.0  # other multipliers
```

Expected at the default settings: baseline 1.85% optimistic (median 0.04 deg,
max 0.24 deg); at 3x the delta-box half-width, 3.85% optimistic (max 0.35 deg).
Optimistic cases average a longer lever arm than the sampled population as a
whole at both box widths; the off-axis-angle difference is pronounced at the
baseline width (0.12 vs 0.07 rad) but nearly disappears at 3x delta (0.08 vs
0.07 rad).

## Reproducing §V.A's sampling-cannot-certify claim

```bash
python experiments/e1_sampling_certification_check.py             # seed 42 (paper text)
python experiments/e1_sampling_certification_check.py --seed 7    # other seeds
```

Phi_minus = max_Theta phi_minus(Theta) is a maximum over the box; Phi_plus is
a minimum. Taking the sample max/min of phi_minus/phi_plus over n draws of
Theta therefore always gives an interval at least as wide as the true robust
cone -- a finite sample can fall short of the true worst-case realization on
each edge, never exceed it. Expected at seed 42: n=5,000 gives sample
[-70.41, 68.10] deg, 18.68 deg wider in aperture than the true cone
([-59.91, 59.91], aperture 119.83 deg); n=20,000 gives 10.10 deg excess. The
exact excess is seed-dependent (it is a single draw, not an average), but the
qualitative point -- the excess shrinks slowly, not proportionally to n -- is
robust across seeds.

## Reproducing E5's rotation-failure horizon

`run_all.py::exp5` reports only pass/fail per trial, not the rotation angle
at which Strategy B (world-frame command) first fails to stick:

```bash
python experiments/e5_rotation_horizon_check.py
```

Must run strategies A, B, C in that order on the same `numpy.random.Generator`
instance (matching `exp5` exactly, which does not reset the RNG between
strategies) -- testing Strategy B in isolation with a fresh seed draws from
the wrong part of the random stream and gives a different, wrong answer.
Expected: A 0.00%, B 96.81% (median rotation-at-failure ~44.1 deg), C 38.42%.

## Reproducing Remark 4 / E4's linear-estimate accuracy claim

`run_all.py::exp4` only evaluates the linear-vs-exact ratio at two box
scales, lambda=0.02 ("small") and lambda=1.0 ("full"). Remark 4 and the E4
paragraph both claim the linear estimate is accurate to within ~5% across
the whole lambda <= 0.1 regime, not just at 0.02:

```bash
python experiments/e4_linear_range_check.py                 # 5 points up to lambda=0.1
python experiments/e4_linear_range_check.py --max-lambda 0.2 --steps 10
```

Expected: relative error grows monotonically from 0.77% (lambda=0.02) to
3.94% (lambda=0.1) -- within the claimed 5% bound at every tested point.

## Validation

`python tests/test_theory.py` checks the theorems against the implementation
(monotonicity in mu/c/alpha including off-axis c, exactness vs a dense delta
grid, the candidate-set bound, the Theorem 3 sticking rotation rate, and
emptiness detection). All 11 tests must pass.

## Notes and caveats

- Experiment E1's `speedup_paper` (225x) is the wall-clock ratio measured once
  on the paper's reference machine against a NumPy-vectorized 21^4 grid
  baseline. `run_all.py` does not re-run that vectorized grid (194,481 points
  x 2,000 geometries would dominate the run time); instead it reports
  `speedup_est`, a quick extrapolation from a scalar per-point evaluation loop.
  `speedup_est` is expected to be much larger than 225x -- it is a cruder,
  deliberately conservative-in-the-other-direction sanity check, not a
  reproduction of the paper's exact figure. Both agree the grid is orders of
  magnitude slower than Algorithm 1; only that qualitative conclusion, plus
  the exact cone edges and the 1.85% vertex-only optimism figure, are load-
  bearing claims in the paper. See the comment above `speedup_est` in
  `experiments/run_all.py::exp1` for the full explanation.
- The paper's "13x slower / 150x slower" figures for the 5^4- and 9^4-point
  polytopic grid baseline (Section V.A) are likewise a one-time wall-clock
  measurement, not reproduced by any script in this package. A manual
  timing of a 5^4/9^4 grid sweep with `cone_angles` against `robust_cone_exact`
  reproduces the same order of magnitude but not identical multipliers --
  same caveat class as `speedup_paper` above: the qualitative conclusion
  (Algorithm 1 is orders of magnitude faster, independent of grid resolution)
  is load-bearing, the exact multipliers are not.
- Experiment E6 in the paper uses PyBullet's rigid-body contact solver, which is
  a *different* contact model from the analytical ellipsoidal limit surface. A
  single PyBullet run instantiates one effective `(mu, c)` realization, so its
  stick-to-slide transition falls between the robust and nominal analytical
  edges. The exact transition angle depends on the engine's friction settings;
  the realized `(mu, c)` are not directly observable from PyBullet alone.
  `pybullet_sweep.py` identifies them explicitly (`identify_effective_mu`, a
  rotation-locked friction probe, and `identify_effective_c`, a bisection
  against `robust_cone_exact`), but this identification is not used as
  validation evidence: the mu-probe has a code-level limitation (a
  `resetBasePositionAndOrientation` side effect that partially masks friction
  sensitivity), and the c-bisection guarantees a near-zero predicted-vs-
  observed residual by construction whenever a solution exists in-bracket,
  independent of whether mu_eff is accurate. See
  `supplementary/pybullet_identification_retraction.md`. The paper relies
  only on the raw observed transition (an unfitted geometric fact) and the
  Simscape mu-sweep (genuinely friction-sensitive, non-circular) for its
  independent-engine friction claims.
- The PyBullet sweep script (`experiments/pybullet_sweep.py`) uses kinematic
  position control (pose reset per step) rather than force/velocity control, and
  classifies stick/slide via a geometric heuristic. These are documented
  limitations; the script serves as a qualitative consistency check, not a
  high-fidelity validation.
- The code was drafted with AI assistance and verified by the authors. Validate
  before use in safety-critical settings.
