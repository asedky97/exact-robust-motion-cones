# Exact Robust Motion Cones for Planar Pushing: Reproducibility Package

[![tests](https://github.com/asedky97/exact-robust-motion-cones/actions/workflows/tests.yml/badge.svg)](https://github.com/asedky97/exact-robust-motion-cones/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

NumPy reference implementation, theory test suite, and experiment drivers for
the IEEE Access submission "Exact Robust Motion Cones for Planar Pushing", by
Ahmed S. M. Elsayed and Sotirios Grammatikos (NTNU). The badge above is live:
it reflects the actual result of the last run of the 11-test theory suite,
the E1–E6 sanity pass, and the five diagnostic scripts below (Remark 1, E1's
median-error/correlation/wider-box, E5's rotation horizon, Remark 4/E4's
λ-range sub-claims, §V.A's sampling-certification claim), executed fresh by
GitHub on a clean machine across 3 operating systems and 3 Python versions on
every push.

## Contents

| Path | Contents |
|------|---------|
| `src/robust_motion_cone.py` | NumPy implementation of Algorithm 1 (exact robust cone) + geometry/uncertainty types |
| `src/pushing_sim.py` | Quasi-static sticking pusher contact model (`resolve_contact`) |
| `experiments/run_all.py` | Experiments E1–E6; writes `experiments/summary.json` |
| `experiments/make_fig7_composite.py` | Assembles Fig. 7 (Simscape independent-engine validation, 2 panels) from the two underlying sweeps' JSON exports; PyBullet is reported in text only, not in this figure; see `REPRODUCIBILITY.md` |
| `experiments/branch_alignment_check.py` | Reproduces Remark 1 (branch-alignment recentre statistics) and the §III.C quadratic-degeneracy note over 10,000 geometries; run with `--wide` for the adversarial 10,000-geometry sweep (full-circle normal angle, lever arms down to 3 cm) also cited in Remark 1; Remark 1 reports on both sweeps combined |
| `experiments/e1_optimism_diagnostics.py` | Reproduces E1's median-error, lever-arm/off-axis correlation, and wider-δ-box sub-claims (§V.A) not covered by `run_all.py` |
| `experiments/e5_rotation_horizon_check.py` | Reproduces E5's ~44° rotation-failure horizon (§V.E) not covered by `run_all.py` |
| `experiments/e4_linear_range_check.py` | Reproduces Remark 4 / E4's "within 5% for λ ≤ 0.1" claim across a range of scales, not just `run_all.py`'s single λ=0.02 point |
| `experiments/e1_sampling_certification_check.py` | Reproduces §V.A's "sampling cannot certify worst-case safety" claim: sample-based max/min estimates of the robust cone edges over n=5,000 and n=20,000 draws, compared against the exact edges |
| `experiments/pybullet_box_sweep.py` (+ `pybullet_box_sweep_world.py` helper) | Reproduces E6(iii): box-spanning independent-solver check; μ swept uniformly across the full box against Zones R/N/S (mirroring E2), rather than one fitted realization |
| `experiments/summary.json` | Canonical headline numbers (reproduced by `run_all.py`) |
| `tests/test_theory.py` | 11 unit tests validating Theorems 1–3 and Proposition 1 |
| `.github/workflows/tests.yml` | Cross-platform CI (Ubuntu/Windows/macOS, Python 3.9/3.11/3.12) |
| `supplementary/validate_matlab.m` | MATLAB reimplementation; cross-checks the cone edges and structural theorems |
| `supplementary/visualize_robust_cone.m` | Animated MATLAB visualization of the robust cone under a parameter sweep (δ, α, μ, c); not part of the validation pipeline, illustrative only |
| `supplementary/` | Simscape Multibody and PyBullet validation scripts + models, Lemma 1 proof, PyBullet identification method and status |
| `REPRODUCIBILITY.md` | Full reproducibility instructions (versions, commands, expected outputs) |

## Quick start

```bash
python -m pip install -r requirements.txt
python -m pytest tests -q          # -> 11 passed
python experiments/run_all.py       # reproduces E1–E6, writes summary.json (~40 s)
```

Expected headline outputs. The geometric/statistical values below are
deterministic or seeded and must match exactly; wall-clock timing numbers
that also appear in `summary.json` (`t_exact_ms`, `t_grid_est_ms`,
`speedup_est`) are hardware- and run-dependent by nature and are omitted
here; they are expected to differ between machines and even between runs
on the same machine:

```
E1 robust=[-59.91, 59.91] nominal=[-85.11, 85.11] optimistic=1.85% speedup 225x
E2 Zone R=0.0%  N=5.24%  S=93.25%
E3 slides=0  spread=9.66deg (bound 33.18deg)  pos_disp=26.4mm (bound 144.8mm)
E4 L_c=155.93/m  L_delta=147.71/m  ratio_full=0.32
E5 A(robust)=0.0%  B(world)=96.81%  C(nominal)=38.42%
E6 transition=±80° (illustrative self-consistency check only; see supplementary/*.m and experiments/pybullet_sweep.py for the paper's actual E6 result)
```

("speedup 225x" is the paper's own reference figure, `speedup_paper`, a
fixed value quoted in the text rather than something `run_all.py`
re-measures; see `REPRODUCIBILITY.md`'s notes on `speedup_est` vs.
`speedup_paper`.)

License: MIT (see `LICENSE`). See `REPRODUCIBILITY.md` for full details.
