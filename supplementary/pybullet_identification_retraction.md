# PyBullet Effective-Parameter Identification: Method and Current Status

**Paper:** *Exact Robust Motion Cones: Guaranteed Sticking Manipulation under Parametric and Geometric Contact Uncertainty.*

This document covers `pybullet_sweep.py`'s attempt to identify each simulated realization's *effective* contact parameters (μ, c) directly, so that a predicted cone can be compared against an observed one for that specific realization: a stronger test, in principle, than checking that an observed transition falls inside the robust/nominal envelope. It records what the approach measures, a code-level issue found in the friction probe and its fix, a separate structural reason the resulting residual can't be read as validation evidence, and what independent-engine evidence from this line of work does stand.

---

## 1. Summary

1. **`identify_effective_mu`'s reported μ_eff = 0.1317 is a measurement-code artifact, not a friction measurement** (Section 3).
2. **Independently of (1), `identify_effective_c`'s bisection guarantees a near-zero predicted-vs-observed residual whenever a solution exists in its search bracket**, regardless of whether μ_eff is correct; the residual is circular by mathematical construction (Section 4).
3. **Unaffected by either issue:** the main PyBullet transition sweep (θ_obs, from `run_push`, a different code path that does not lock rotation), all Python theory/experiments (E1–E5), the Simscape validation (continuous penalty contact, not kinematic teleport-and-reset: a structurally different contact treatment from the bug class in Section 3), and the independent MATLAB reimplementation.

The paper accordingly does not report a fitted-parameter PyBullet match. Its independent-engine friction validation rests on the Simscape sweeps (`supplementary/simscape_validation_status.md`, Labs 1–2, both non-circular); the raw PyBullet transition angle (θ_obs) stands only as a qualitative geometric observation, not a parameter-identified quantitative test.

---

## 2. Method: identifying (μ, c) per realization

`pybullet_sweep.py` identifies each realization's effective contact parameters with two probes, run as separate experiments from the main motion-cone sweep so identification does not beg the question being tested:

**Effective friction, μ_eff (`identify_effective_mu`).** Push through the slider's centre of mass, so the contact force produces zero torque and the body translates without rotating. In this regime the stick/slide boundary is the friction cone itself: μ_eff = tan(θ_slide), where θ_slide is the push angle (from the contact normal) at which lateral slip onsets.

**Effective contact arm, c_eff (`identify_effective_c`).** With μ_eff fixed, push at the actual contact point along the normal. The off-centre force produces both translation and rotation; the ratio of steady-state angular to translational velocity is set by the limit-surface coupling, and inverting it (via `robust_motion_cone.py`) recovers the effective contact arm.

**Predicted cone.** With (μ_eff, c_eff) identified, the predicted half-angle for the realization is the *nominal* cone at those parameters (a zero-width box), evaluated with `robust_cone_exact`.

**Observed cone.** The ordinary motion-cone sweep at the same geometry gives the observed transition θ_obs: the largest push angle with sustained sticking.

The deliverable is the residual |θ_obs − θ_pred|: a direct test of the cone *equation*, rather than only checking θ_obs ∈ [robust, nominal].

---

## 3. Why μ_eff is not a valid friction measurement

`identify_effective_mu` calls a rotation-locked push probe that re-centers the board's orientation to identity every simulation step:

```python
xpos, _ = p.getBasePositionAndOrientation(board)
p.resetBasePositionAndOrientation(board, [xpos[0], xpos[1], xpos[2]], [0, 0, 0, 1])
```

`resetBasePositionAndOrientation` zeroes **both linear and angular velocity** as a side effect, even when position is passed back unchanged and only orientation is being corrected:

```
velocity before reset: ((1.0, 2.0, 3.0), (0.1, 0.2, 0.3))
velocity AFTER resetBasePositionAndOrientation (same pos, identity orn): ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
```

Because this call runs every step immediately before the relative tangential velocity is sampled, `getBaseVelocity(board)` always reads zero translation for the board, so the probe's `vrel_t = pusher_v − v_com` collapses to `pusher_v` itself, independent of the board, independent of contact, independent of μ.

**Consequence.** Sweeping pusher friction from 0.01 to 5.0 (500×) and board friction from 0.01 to 0.9, the probe's reported slide-onset angle does not change at all; it always lands at the angle where `sin(θ) = 0.15` (the classification tolerance) crosses, i.e. θ = arcsin(0.15) = 8.627°, floor-rounded to the nearest tested grid point. At a 2.5°-step grid, that point is 7.5°, and tan(7.5°) = 0.1317: exactly the value reported as μ_eff, for any underlying friction. The number was never measuring friction.

**Fix.** Capture linear velocity immediately before the orientation-lock reset and restore it afterward (angular velocity is left zeroed, consistent with the intent of locking rotation):

```python
xpos, _ = p.getBasePositionAndOrientation(board)
vlin_pre, _ = p.getBaseVelocity(board)
p.resetBasePositionAndOrientation(board, [xpos[0], xpos[1], xpos[2]], [0, 0, 0, 1])
p.resetBaseVelocity(board, linearVelocity=vlin_pre, angularVelocity=[0, 0, 0])
```

**Effect of the fix** (fresh-world tests, pusher friction swept against the fixed probe):

| pusher μ (set) | θ_slide (fixed probe) | μ_eff = tan(θ_slide) | theory atan(μ) |
|---|---|---|---|
| 0.05 | 11.0° | 0.194 | 2.9° |
| 0.25 | 19.0° | 0.344 | 14.0° |
| 0.45 | 19.0° | 0.344 | 24.2° |
| 0.65 | 19.0° | 0.344 | 33.0° |
| 1.50 | 19.0° | 0.344 | 56.3° |

The fix restores partial friction sensitivity (μ = 0.05 now differs from higher values, as physics demands) but the result saturates for μ ≥ 0.25. Achieved contact normal force at μ = 0.25 vs μ = 1.5 (mean 3.60 N vs 4.05 N over the sampling window) shows Fn is not strictly capped, so normal-force saturation alone does not fully explain the plateau; some additional interaction between the kinematic teleportation step size, contact stiffness, and the fixed 15%-of-push-speed classification tolerance is likely also involved. **This residual discrepancy is not fully diagnosed.** The fix (a real correctness improvement, kept in `pybullet_sweep.py`) is not treated as producing a validated friction measurement, and a code comment to that effect sits at the call site, rather than presenting a still-imperfect number as fully validated or spending unbounded effort chasing a fully friction-linear kinematic-pusher probe when a genuinely friction-sensitive independent-engine result (the Simscape μ-sweep) already exists.

---

## 4. Why the c_eff residual is circular by construction

`identify_effective_c` does not call the rotation-locked probe above at all, so it is unaffected by Section 3. It is a pure bisection:

```python
f_lo = _predicted_edge(geom, mu_eff, c_lo) - theta_obs_rad
f_hi = _predicted_edge(geom, mu_eff, c_hi) - theta_obs_rad
# ... bisect until _predicted_edge(geom, mu_eff, c_mid) == theta_obs_rad
```

For any fixed μ_eff (correct or not) and any observed transition that admits a solution c ∈ [c_lo, c_hi], this converges to a c_eff such that the predicted edge matches the observation to within the bisection tolerance, by construction of the search, not by empirical agreement. A near-zero predicted-vs-observed residual is therefore guaranteed whenever a solution exists in-bracket; it is not evidence that the ellipsoidal model matches independent physics at that point. This holds regardless of Section 3's bug and is the more fundamental reason the "predicted = observed" framing is not used as a headline validation claim.

---

## 5. What stands as independent-engine evidence

| Test | Fits a free parameter to match an observation? | Status |
|---|---|---|
| Simscape 44-angle sweep at fixed nominal (μ₀, c₀) = (0.45, 0.06) | No | Solid: 84.1% (37/44); see `simscape_validation_status.md` |
| Simscape μ-sweep, μ ∈ {0.20, ..., 0.60} | No | Solid: exact atan(μ) recovery (0.00° error at every point) |
| MATLAB independent reimplementation | No (separate codebase, same closed-form scenario) | Solid: all assertions pass |
| PyBullet θ_obs (main sweep, `run_push`) | No | Stands as a raw geometric observation; not compared against a fitted prediction |
| PyBullet μ_eff (`identify_effective_mu`) | N/A (measurement, not a fit) | Not used as a validated friction number; see Section 3 |
| PyBullet c_eff + predicted-vs-observed residual | Yes, by bisection to the observation | Not used as validation evidence (see Section 4); the underlying geometric fact θ_obs < θ_nominal remains true, just not usable as a quantitative model check |

The paper's independent-engine claims rest on the first three rows.
