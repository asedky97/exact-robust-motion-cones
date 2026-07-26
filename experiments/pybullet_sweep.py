#!/usr/bin/env python3
"""
Independent-engine validation (paper Experiment E6).

Runs a dense push-angle sweep in PyBullet, an independent rigid-body physics
engine that resolves stick/slide through its own contact solver rather than the
ellipsoidal limit surface used by the analytical theory. This is an OUT-OF-MODEL
check: agreement confirms the analytical pushing map captures the contact
physics, and that the realized cone for any single parameter realization lies
strictly inside the nominal boundary.

PyBullet is NOT required by the core algorithm or the other experiments; it is
only used here. Install with:  pip install pybullet

Usage:
    python pybullet_sweep.py                 # 45-angle sweep, headless
    python pybullet_sweep.py --gui           # show the simulation window
    python pybullet_sweep.py --angles 90     # finer sweep

Output:
    ../figures/fig_pybullet_sweep.png/.pdf
    ../experiments/pybullet_results.json     (per-angle stick/slide + efficiency)

LIMITATIONS:
- The pusher is controlled by resetting its pose each simulation step
  (kinematic position control) rather than by a velocity servo or force
  control. This may excite transient dynamics at each step. A more faithful
  implementation would use a prismatic joint with `setJointMotorControl2`.
- The stick/slide classifier uses a geometric heuristic (translation along
  push > 50% of distance AND lateral slip < 20% of along). This heuristic
  is reasonable for quasi-static planar pushing but is not equivalent to the
  analytical friction-cone condition used in the paper. The results should
  be interpreted as a consistency check, not a direct validation.
- The realized effective (mu, c) are identified per realization by two short
  probes (Q6): a rotation-locked friction probe for mu_eff = tan(theta_slide),
  and a 1-D bisection on the contact arm for c_eff. The observed transition is
  then compared to the cone predicted by robust_cone_exact at (mu_eff, c_eff),
  giving a direct predicted-vs-observed test (residual reported to stdout and
  pybullet_results.json). The stick/slide classifier remains a geometric
  heuristic; the identification converts the envelope check into a direct
  cone-equation test but does not make the classifier first-principles.
- Bullet's kinematic (teleported) pusher typically realizes an effective
  friction mu_eff BELOW the paper's uncertainty box (mu in [0.25, 0.65]); the
  tangential drag of a reset-per-step contact is under-resolved. The observed
  cone is then correctly NARROWER than the robust envelope (a sub-box friction
  gives a sub-robust cone). The Q6 deliverable is therefore the predicted-vs-
  observed match at the *realization's identified* (mu_eff, c_eff) -- not
  envelope membership -- and the Simscape mu-sweep independently confirms the
  cone equals atan(mu) across mu in {0.20..0.60}, spanning the box.
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from robust_motion_cone import ContactGeometry, UncertaintyBox, cone_angles, robust_cone_exact  # noqa

# Paper scenario (board + table); tune friction to a representative realization.
BOARD_HALF = (0.30, 0.05, 0.0125)      # 0.60 x 0.10 x 0.025 m
PUSH_SPEED = 0.003                      # 3 mm/s (quasi-static)
PUSH_DIST = 0.040                       # 40 mm
TABLE_FRICTION = 0.30                   # representative effective friction (board-table)
BOARD_MASS = 1.2
PUSHER_FRICTION = 0.45                   # pusher-board contact friction (paper mu0); set
                                        # explicitly so Q6 identification has a known nominal.


def build_world(gui=False):
    import pybullet as p
    import pybullet_data
    cid = p.connect(p.GUI if gui else p.DIRECT)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.setGravity(0, 0, -9.81)
    p.loadURDF("plane.urdf")
    # board as a box
    col = p.createCollisionShape(p.GEOM_BOX, halfExtents=BOARD_HALF)
    vis = p.createVisualShape(p.GEOM_BOX, halfExtents=BOARD_HALF,
                              rgbaColor=[0.76, 0.64, 0.42, 1])
    board = p.createMultiBody(BOARD_MASS, col, vis,
                              basePosition=[0, 0, BOARD_HALF[2] + 1e-3])
    p.changeDynamics(board, -1, lateralFriction=TABLE_FRICTION)
    # pusher: a small cylinder kinematically controlled (teleported each step
    # at the push velocity).  A kinematic pusher develops a reliable, repeatable
    # normal contact; the vrel_t contact-point classifier (in run_push) measures
    # the true stick/slide boundary so the identified (mu_eff, c_eff) and the
    # observed cone are self-consistent regardless of the realized friction.
    pcol = p.createCollisionShape(p.GEOM_CYLINDER, radius=0.005, height=0.03)
    pusher = p.createMultiBody(0, pcol, basePosition=[-0.30, 0, 0.015])
    p.changeDynamics(pusher, -1, lateralFriction=PUSHER_FRICTION)
    return cid, board, pusher


def run_push(board, pusher, angle_deg, steps=2400):
    """Quasi-static push: pusher advances PUSH_DIST over `steps` steps.

    Returns (mode, eff) where mode is 'stick'/'slide' and eff is the translation
    efficiency [%] along the push direction. The pusher is kinematically swept
    a total distance PUSH_DIST (40 mm) over the run so the stick heuristic
    (along > 0.5*PUSH_DIST, lateral < 20% along) is meaningful.
    """
    import pybullet as p
    # reset board pose
    p.resetBasePositionAndOrientation(board, [0, 0, BOARD_HALF[2] + 1e-3], [0, 0, 0, 1])
    p.resetBaseVelocity(board, [0, 0, 0], [0, 0, 0])
    start = np.array([-0.30, 0.0, 0.015])
    p.resetBasePositionAndOrientation(pusher, start, [0, 0, 0, 1])
    a = np.deg2rad(angle_deg)
    direction = np.array([np.cos(a), np.sin(a), 0.0])
    dt = 1.0 / 240
    p.setTimeStep(dt)
    pos0, _ = p.getBasePositionAndOrientation(board)
    per_step = PUSH_DIST / steps   # pusher displacement per step [m]
    push_speed = per_step / dt     # m/s; set the kinematic pusher's linear velocity so
                                  # Bullet's friction solver sees the relative tangential
                                  # velocity and develops a proper tangential drag (without
                                  # this a teleported kinematic body has ~zero velocity and
                                  # the board is only shoved along the normal).
    p0_body = np.array([-0.30, 0.0, 0.0])   # contact point on -x face, body frame (COM at origin)
    n_ss = max(1, steps // 5)               # last 20% of steps -> steady state
    vrel_t_samples = []
    for i in range(1, steps + 1):
        new = start + direction * (per_step * i)
        p.resetBasePositionAndOrientation(pusher, new.tolist(), [0, 0, 0, 1])
        p.resetBaseVelocity(pusher, [direction[0]*push_speed, direction[1]*push_speed, 0.0], [0,0,0])
        p.stepSimulation()
        if i > steps - n_ss:
            vlin, vang = p.getBaseVelocity(board)
            _, orn = p.getBasePositionAndOrientation(board)
            yaw = p.getEulerFromQuaternion(orn)[2]
            v_com = np.array(vlin[:2])
            wz = vang[2]
            cy, sy = np.cos(yaw), np.sin(yaw)
            r_world = np.array([cy*p0_body[0] - sy*p0_body[1], sy*p0_body[0] + cy*p0_body[1]])
            v_cp = v_com + np.array([-wz*r_world[1], wz*r_world[0]])   # v_com + omega x r (2D)
            t_world = np.array([-sy, cy])                            # body tangent [0,1] -> world
            pusher_v = np.array([direction[0]*push_speed, direction[1]*push_speed])
            vrel_t_samples.append(float(np.dot(pusher_v - v_cp, t_world)))
    pos1, orn1 = p.getBasePositionAndOrientation(board)
    disp = np.array(pos1[:2]) - np.array(pos0[:2])
    along = float(np.dot(disp, direction[:2]))
    eff = 100 * max(0.0, along) / PUSH_DIST
    # Motion-cone metric: relative tangential velocity at the contact point.
    # Sticking (pusher velocity inside the cone) => |vrel_t| ~ 0; sliding (outside)
    # => |vrel_t| > 0.  A rotating board that STICKS moves its contact point with the
    # pusher even though its COM drifts off the push line, so the old COM-translation
    # heuristic mis-labelled such cases as 'slide'.  vrel_t measures the contact directly.
    vrel_t_ss = float(np.mean(np.abs(vrel_t_samples))) if vrel_t_samples else 0.0
    tol = 0.15 * push_speed
    mode = "stick" if vrel_t_ss < tol else "slide"
    return mode, eff


def _run_push_probe(board, pusher, angle_deg, contact_y=0.0, lock_rotation=False, steps=1800):
    """Push variant for contact-parameter identification (Q6).

    contact_y : y-coordinate of the contact point on the -x face; changing it
                changes the lever arm about the COM and hence the rotation
                demand (used to isolate the contact arm).
    lock_rotation : if True the board orientation is reset to identity every
                step, enforcing PURE TRANSLATION. With rotation suppressed the
                stick/slide onset is governed by the contact friction alone,
                so mu_eff = tan(theta_slide). This is a SEPARATE experiment from
                the motion-cone sweep (no torque -> no dependence on c), so the
                identification does not beg the question.
    """
    import pybullet as p
    p.resetBasePositionAndOrientation(board, [0, 0, BOARD_HALF[2] + 1e-3], [0, 0, 0, 1])
    p.resetBaseVelocity(board, [0, 0, 0], [0, 0, 0])
    start = np.array([-0.30, contact_y, 0.015])
    p.resetBasePositionAndOrientation(pusher, start, [0, 0, 0, 1])
    a = np.deg2rad(angle_deg)
    direction = np.array([np.cos(a), np.sin(a), 0.0])
    dt = 1.0 / 240
    p.setTimeStep(dt)
    pos0, _ = p.getBasePositionAndOrientation(board)
    per_step = PUSH_DIST / steps
    push_speed = per_step / dt     # see run_push: lets the friction solver see the pusher velocity
    p0_body = np.array([-0.30, contact_y, 0.0])   # contact point on -x face, body frame
    n_ss = max(1, steps // 5)
    vrel_t_samples = []
    for i in range(1, steps + 1):
        new = start + direction * (per_step * i)
        p.resetBasePositionAndOrientation(pusher, new.tolist(), [0, 0, 0, 1])
        p.resetBaseVelocity(pusher, [direction[0]*push_speed, direction[1]*push_speed, 0.0], [0,0,0])
        p.stepSimulation()
        if lock_rotation:
            # NOTE: resetBasePositionAndOrientation silently zeroes both
            # linear and angular velocity as a side effect, even
            # when position/orientation are otherwise unchanged (verified via
            # getBaseVelocity before/after on an isolated test body). Calling
            # it every step to suppress rotation was therefore also erasing
            # the board's translational velocity before it could be read
            # below, making vrel_t_ss collapse to the commanded pusher
            # velocity itself -- independent of the contact friction -- and
            # mu_eff = tan(theta_slide) into a pure classification-threshold
            # artifact (tan(arcsin(0.15)) rounded to the nearest tested grid
            # point) rather than a friction measurement. Fix: capture linear
            # velocity before the reset and restore it after; zero angular
            # velocity explicitly since orientation is being forced to 0 anyway.
            xpos, _ = p.getBasePositionAndOrientation(board)
            vlin_pre, _ = p.getBaseVelocity(board)
            p.resetBasePositionAndOrientation(board, [xpos[0], xpos[1], xpos[2]], [0, 0, 0, 1])
            p.resetBaseVelocity(board, linearVelocity=vlin_pre, angularVelocity=[0, 0, 0])
        if i > steps - n_ss:
            vlin, vang = p.getBaseVelocity(board)
            _, orn = p.getBasePositionAndOrientation(board)
            yaw = p.getEulerFromQuaternion(orn)[2]
            v_com = np.array(vlin[:2]); wz = vang[2]
            cy, sy = np.cos(yaw), np.sin(yaw)
            r_world = np.array([cy*p0_body[0] - sy*p0_body[1], sy*p0_body[0] + cy*p0_body[1]])
            v_cp = v_com + np.array([-wz*r_world[1], wz*r_world[0]])
            t_world = np.array([-sy, cy])
            pusher_v = np.array([direction[0]*push_speed, direction[1]*push_speed])
            vrel_t_samples.append(float(np.dot(pusher_v - v_cp, t_world)))
    pos1, _ = p.getBasePositionAndOrientation(board)
    disp = np.array(pos1[:2]) - np.array(pos0[:2])
    along = float(np.dot(disp, direction[:2]))
    lateral = float(np.linalg.norm(disp - along * direction[:2]))
    vrel_t_ss = float(np.mean(np.abs(vrel_t_samples))) if vrel_t_samples else 0.0
    tol = 0.15 * push_speed
    mode = "stick" if vrel_t_ss < tol else "slide"
    return mode, along, lateral


def identify_effective_mu(board, pusher, verbose=True):
    """Identify the effective pusher-board friction mu_eff (Q6).

    Rotation-locked friction probe: push at the face centre with the board's
    rotation kinematically suppressed, sweep the push angle, and take
        mu_eff = tan(theta_slide)
    where theta_slide is the largest angle that still sticks. With no rotation
    the stick/slide boundary is the friction cone, independent of the contact
    arm c, so this isolates mu.

    CAUTION (see supplementary/pybullet_identification_retraction.md): a
    resetBasePositionAndOrientation velocity-clobbering bug in the
    lock_rotation path (fixed above) previously made this
    function's output completely independent of the true friction
    coefficient. The fix restores partial sensitivity but the result still
    saturates for mu >= ~0.25 for reasons not fully diagnosed (achieved
    contact normal force differs only ~12% between mu=0.25 and mu=1.5, so
    normal-force saturation alone does not fully explain it). Do not treat
    this function's output as a validated friction measurement without
    further investigation; the paper's independent-engine friction
    validation relies on the Simscape mu-sweep instead (genuinely
    friction-sensitive, exact atan(mu) recovery, no fitting).
    """
    angles = np.arange(5, 46, 2.5)            # 5 .. 42.5 deg
    stick_angles = []
    for ad in angles:
        mode, _, _ = _run_push_probe(board, pusher, ad, contact_y=0.0, lock_rotation=True)
        if mode == "stick":
            stick_angles.append(ad)
        elif verbose:
            print(f"  mu-probe {ad:5.1f} deg -> {mode}")
    if not stick_angles:
        return None
    theta_slide = max(stick_angles)
    mu_eff = np.tan(np.deg2rad(theta_slide))
    if verbose:
        print(f"  mu-probe: theta_slide = {theta_slide:.1f} deg -> mu_eff = {mu_eff:.3f}")
    return float(mu_eff)


def _predicted_edge(geom, mu, c):
    """Predicted nominal cone half-angle [rad] for identified (mu, c), zero box."""
    zb = UncertaintyBox(mu=(mu, mu), c=(c, c), alpha=(0.0, 0.0), delta=(0.0, 0.0))
    _, hi, _, _ = robust_cone_exact(geom, zb)
    return hi


def identify_effective_c(geom, mu_eff, theta_obs_rad, c_lo=0.03, c_hi=0.12, iters=40, verbose=True):
    """Identify the effective contact arm c_eff (Q6).

    With mu_eff known, the predicted motion-cone half-angle is monotone in c
    (Theorem 1(ii)). We bisection-search c so that the predicted edge matches the
    observed transition theta_obs:
        robust_cone_exact(geom, zero-box@(mu_eff, c_eff, 0, 0)).hi == theta_obs.
    """
    f_lo = _predicted_edge(geom, mu_eff, c_lo) - theta_obs_rad
    f_hi = _predicted_edge(geom, mu_eff, c_hi) - theta_obs_rad
    if f_lo * f_hi > 0:
        if verbose:
            print(f"  c-probe: no sign change on [{c_lo},{c_hi}] "
                  f"(edges {np.rad2deg(f_lo+theta_obs_rad):.1f},{np.rad2deg(f_hi+theta_obs_rad):.1f} vs obs {np.rad2deg(theta_obs_rad):.1f})")
        return None
    for _ in range(iters):
        c_mid = 0.5 * (c_lo + c_hi)
        f_mid = _predicted_edge(geom, mu_eff, c_mid) - theta_obs_rad
        if f_lo * f_mid <= 0:
            c_hi, f_hi = c_mid, f_mid
        else:
            c_lo, f_lo = c_mid, f_mid
    c_eff = 0.5 * (c_lo + c_hi)
    if verbose:
        print(f"  c-probe: c_eff = {c_eff:.4f} m "
              f"(predicted edge {np.rad2deg(_predicted_edge(geom, mu_eff, c_eff)):.2f} deg)")
    return float(c_eff)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gui", action="store_true")
    ap.add_argument("--angles", type=int, default=45)
    args = ap.parse_args()

    try:
        import pybullet  # noqa
    except ImportError:
        sys.exit("PyBullet not installed. Run: pip install pybullet")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    box = UncertaintyBox(mu=(0.25, 0.65), c=(0.045, 0.075),
                         alpha=(np.deg2rad(-6), np.deg2rad(6)), delta=(-0.015, 0.015))
    geom = ContactGeometry(p0=np.array([-0.30, 0.0]), n0=np.array([1.0, 0.0]))
    lo, hi, ref, _ = robust_cone_exact(geom, box)
    mu0, c0, a0, d0 = box.nominal()
    lon, hin = cone_angles(geom, mu0, c0, a0, d0, ref)
    ROB, NOM = np.rad2deg(hi), np.rad2deg(hin)

    cid, board, pusher = build_world(args.gui)
    angles = np.linspace(-88, 88, args.angles)
    modes, effs = [], []
    t0 = time.perf_counter()
    for ad in angles:
        m, e = run_push(board, pusher, ad)
        modes.append(m); effs.append(e)
    print(f"swept {len(angles)} angles in {time.perf_counter()-t0:.1f}s")

    stick = [a for a, m in zip(angles, modes) if m == "stick"]
    transition = max(stick) if stick else None

    # ---- Q6: effective (mu, c) identification + predicted-vs-observed ----
    mu_eff = identify_effective_mu(board, pusher, verbose=True)
    c_eff = None
    theta_pred_deg = None
    residual_deg = None
    if transition is not None and mu_eff is not None:
        c_eff = identify_effective_c(geom, mu_eff, np.deg2rad(transition), verbose=True)
        if c_eff is not None:
            theta_pred_deg = float(np.rad2deg(_predicted_edge(geom, mu_eff, c_eff)))
            residual_deg = float(abs(transition - theta_pred_deg))
    print("\n--- Q6 predicted-vs-observed ---")
    print(f"  mu_eff            = {mu_eff:.4f}  (nominal mu0 = {mu0:.2f})" if mu_eff else "  mu_eff = <not identified>")
    print(f"  c_eff             = {c_eff:.4f} m (nominal c0 = {c0:.3f})" if c_eff else "  c_eff  = <not identified>")
    print(f"  observed  theta   = +/-{transition:.2f} deg" if transition else "  observed  theta   = <none>")
    print(f"  predicted theta   = +/-{theta_pred_deg:.2f} deg" if theta_pred_deg else "  predicted theta   = <none>")
    print(f"  |obs - pred|      = {residual_deg:.2f} deg" if residual_deg is not None else "  residual = <n/a>")
    print(f"  envelope [robust {ROB:.2f}, nominal {NOM:.2f}] deg "
          f"-> obs {'inside' if (transition and ROB <= transition <= NOM) else '?'} envelope")

    res = {"transition_deg": transition, "robust_edge_deg": round(ROB, 2),
           "nominal_edge_deg": round(NOM, 2),
           "mu_eff": mu_eff, "c_eff": c_eff,
           "predicted_edge_deg": theta_pred_deg, "residual_deg": residual_deg,
           "nominal_mu0": mu0, "nominal_c0": c0,
           "angles": angles.tolist(), "modes": modes, "efficiency": effs}
    (ROOT / "experiments" / "pybullet_results.json").write_text(json.dumps(res, indent=2))

    # Two-panel figure: (a) efficiency sweep with all edges incl. predicted;
    # (b) predicted-vs-observed cone half-angle bar comparison.
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9, 3.2))
    cols = ["#2a8c4a" if m == "stick" else "#c0392b" for m in modes]
    ax.scatter(angles, effs, c=cols, s=20)
    for v, c in [(ROB, "#b8002e"), (-ROB, "#b8002e"), (NOM, "#5b8bd0"), (-NOM, "#5b8bd0")]:
        ax.axvline(v, color=c, lw=1, ls="-" if c == "#b8002e" else "--")
    if theta_pred_deg is not None:
        ax.axvline(theta_pred_deg, color="#000000", lw=1.4, ls=":")
        ax.axvline(-theta_pred_deg, color="#000000", lw=1.4, ls=":")
    ax.set_xlabel("push angle (deg)"); ax.set_ylabel("translation efficiency (%)")
    ax.set_title(f"PyBullet sweep: obs transition +/-{transition:.0f} deg"
                 if transition else "PyBullet sweep")

    labels = ["robust\n(+/-)", "predicted\n(Q6)", "observed", "nominal\n(+/-)"]
    vals = [ROB, theta_pred_deg if theta_pred_deg else 0,
            transition if transition else 0, NOM]
    barcols = ["#b8002e", "#000000", "#2a8c4a", "#5b8bd0"]
    ax2.bar(labels, vals, color=barcols)
    for xi, v in zip(range(len(vals)), vals):
        if v:
            ax2.text(xi, v + 1, f"{v:.1f}", ha="center", fontsize=8)
    ax2.set_ylabel("cone half-angle (deg)")
    ax2.set_title("Predicted vs observed (Q6)")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "fig_pybullet_sweep.png", dpi=300)
    fig.savefig(ROOT / "figures" / "fig_pybullet_sweep.pdf")
    print(f"transition +/-{transition} deg  (robust +/-{ROB:.1f}, nominal +/-{NOM:.1f})"
          if transition else "transition <none>")


if __name__ == "__main__":
    main()
