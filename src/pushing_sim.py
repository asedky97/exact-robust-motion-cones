"""
Quasi-static pusher-slider simulator (ellipsoidal limit surface).

Given a commanded pusher velocity u_p at body-frame contact point p with inward
normal n, friction mu, and limit-surface length c, resolves the contact mode:

  Sticking:  f_dir = M^{-1} u_p lies inside the friction cone
             => twist v = (f_x, f_y, s*tau(f)) with contact velocity exactly u_p.
  Sliding:   contact force saturates on a friction-cone edge f_e; twist scaled so
             the normal velocity component matches the pusher (contact maintained).
  Separation: u_p . n <= 0 (pulling away).

Used by experiments.py for Monte Carlo validation. All quantities body-frame
unless noted; world pose integrated externally.
"""

from __future__ import annotations
import numpy as np
from robust_motion_cone import perp, cross2, friction_edges, push_map

def resolve_contact(p, n, mu, c, u_p):
    """Return (mode, twist) for commanded contact velocity u_p.
    mode in {'stick', 'slide+', 'slide-', 'separate'};
    twist = (vx, vy, omega) body frame."""
    p, n, u_p = map(np.asarray, (p, n, u_p))
    s = 1.0 / c**2
    M = push_map(p, c)
    f_dir = np.linalg.solve(M, u_p)
    fm, fp = friction_edges(n, mu)
    # Sticking iff the required force lies inside the friction cone (CCW sector
    # from fm to fp) and is compressive. NOTE: separation is a FORCE condition
    # (tensile normal force required), never a condition on the velocity
    # direction -- under sticking the contact point may move in any direction.
    inside = (cross2(fm, f_dir) >= -1e-12 and cross2(f_dir, fp) >= -1e-12
              and float(f_dir @ n) > 0)
    if inside:
        tau = cross2(p, f_dir)
        return "stick", np.array([f_dir[0], f_dir[1], s * tau])
    # Sliding: contact force saturates on the friction-cone edge nearer f_dir;
    # twist scaled so normal velocity components match (contact maintained).
    f_e, mode = (fp, "slide+") if cross2(f_dir, fp) < 0 else (fm, "slide-")
    u_e = M @ f_e
    denom = float(u_e @ n)
    kappa = float(u_p @ n) / denom if abs(denom) > 1e-12 else -1.0
    if kappa <= 0:
        return "separate", np.zeros(3)
    tau = cross2(p, f_e)
    return mode, kappa * np.array([f_e[0], f_e[1], s * tau])


def integrate_push(pose, p, n, mu, c, u_p_world, dt, steps):
    """Integrate a constant world-frame pusher velocity for `steps` of size dt.
    pose = (x, y, psi) world pose of the object CoM. The pusher is assumed to
    track the (moving) contact point; contact location p, normal n fixed in
    the body frame (sticking-intent commands). Returns trajectory and modes."""
    traj, modes = [np.array(pose, dtype=float)], []
    pose = np.array(pose, dtype=float)
    for _ in range(steps):
        R = np.array([[np.cos(pose[2]), -np.sin(pose[2])],
                      [np.sin(pose[2]),  np.cos(pose[2])]])
        u_body = R.T @ np.asarray(u_p_world)
        mode, v = resolve_contact(p, n, mu, c, u_body)
        modes.append(mode)
        if mode == "separate":
            break
        # body twist -> world pose rate
        pose = pose + dt * np.array([*(R @ v[:2]), v[2]])
        traj.append(pose.copy())
    return np.array(traj), modes
