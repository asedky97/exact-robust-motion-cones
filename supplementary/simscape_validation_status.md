# Simscape Multibody Validation

**Files:** `supplementary/robust_motion_cone_simscape.m` (motion-cone classification sweep), `supplementary/cone_simscape_sweep.m` (friction-cone sweep), `supplementary/push_sim.slx` (model, built programmatically by the main script)

**Environment:** MATLAB R2024b with Simscape Multibody.

---

## Model

The model is built programmatically (no GUI) via `add_block`/`add_line`:

- **Pusher body**: Brick geometry (0.06 × 0.04 × 0.02 m), actuated by a Prismatic Joint driven by a prescribed velocity (set per angle via `set_param`).
- **Slider body**: Brick geometry (0.5 × 0.3 × 0.02 m), on a Planar Joint with base damping (X: 20 N·s/m, Y: 20 N·s/m, Rz: 0.5 N·m·s/rad) and a 6-DOF sensor for position, velocity, and rotation.
- **Contact offset**: Rigid Transform shifts the pusher frame to the slider's left face (`contact_arm_x = -0.25 m`).
- **External Force and Torque (EFT)** applies contact forces at the offset point (Fx normal, Fy tangential, Tz torque).
- **Contact force law** (continuous penalty model, tuned to hold steady-state normal force through the push):
  ```
  GapX = Px - Sx,  vrelX = Vx_pusher - Vx_slider
  Fn   = max(Fmin, Kx·(GapX+pre_pen) + Dx·vrelX)   : pre-compression + F_min floor
  GapY = Py - Sy,  vrelY = Vy_pusher - Vy_slider
  Ft   = sat(Kt·GapY + Dt·vrelY, ±μ·Fn)             : regularized Coulomb
  Tau  = contact_arm_x · Ft                          : torque from edge contact
  ```
  The pre-compression offset (`pre_pen = 0.2 mm`) and minimum normal force floor (`Fmin = 5 N`) keep the contact engaged through the full push rather than letting the penalty spring relax to zero at steady state; with μ = 0.45 the sustained friction limit μ·Fmin = 2.25 N exceeds the base damping drag, so sticking is physically sustainable inside the cone.
- **Parameters:** Kx = Kt = 10000 N/m, Dx = Dt = 100 N·s/m, `pre_pen = 0.2 mm`, `Fmin = 5 N`.
- **Simulation:** ode45, 1.0 s, MaxStep 0.005, ZeroCross off.
- **Classifier:** stick/slide is decided from the contact-point relative tangential velocity `vrel_t = vy − (vy_ss − contact_arm·drz_ss)`: the physically correct metric, since the slider's translation direction alone is confounded by its rotation. Sticking ⟹ |vrel_t| ≈ 0; sliding ⟹ |vrel_t| > 0. The tolerance is adaptive: `stick_tol = max(5°, ½·g_eff)`, where `g_eff` is the effective friction-cone half-angle read from the simulation's own direction saturation at large push angles, rather than a fixed threshold.

## Lab 1: motion-cone classification (`robust_motion_cone_simscape.m`)

Sweeps 44 angles (−88° to −4° and +4° to +88° in 4° steps, skipping 0°) at μ = 0.45, comparing the simulated stick/slide outcome against the analytical robust-cone prediction (nominal [−85.11°, +85.11°], robust [−59.91°, +59.91°]).

```
Accuracy:  84.1% (37/44)
Precision: 88.9%
Recall:    85.7%
```

Running the script produces `figures/simscape_validation.png` (direction with within-run error bars, speed, speed-ratio, classification, Rz rotation, Fn normal-force diagnostic); not tracked in the repo.

## Lab 2: friction-cone recovery (`cone_simscape_sweep.m`)

A single reusable model is rebuilt per friction value (a unique model name per μ forces a fresh Simulink compile so the block's gain actually updates) and swept over μ ∈ {0.20, 0.30, 0.40, 0.50, 0.60} with a clean one-sided spring-damper contact (no pre-compression or force floor: the pusher's own penetration sets the normal force, so the friction cone emerges without an external offset biasing it).

| μ | theory atan(μ) | sim | err |
|----|----------------|-----|-----|
| 0.20 | 11.31° | 11.31° | 0.00° |
| 0.30 | 16.70° | 16.70° | 0.00° |
| 0.40 | 21.80° | 21.80° | 0.00° |
| 0.50 | 26.57° | 26.57° | 0.00° |
| 0.60 | 30.96° | 30.96° | 0.00° |

Running the script produces `figures/simscape_cone_halfangle_vs_mu.png`; not tracked in the repo.

## What this does and doesn't establish

Labs 1 and 2 are both non-circular: neither fits a free parameter to match an observation. Lab 1 checks classification accuracy against a fixed analytical prediction; Lab 2 recovers a known closed-form relationship (the friction-cone half-angle) exactly, across a range of friction values, from an independent contact solver (Simscape's regularized-penalty law, a different numerical treatment from the analytical ellipsoidal limit surface). Together they are the paper's independent-engine friction validation.

A separate, PyBullet-based attempt to identify each realization's *effective* (μ, c) and compare a predicted cone directly against an observed one is not used as quantitative validation evidence; see `supplementary/pybullet_identification_retraction.md` for why, and for what from that PyBullet work does still stand.
