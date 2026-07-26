"""
World/push helpers for pybullet_box_sweep.py, factored out of pybullet_sweep.py
so mu can be set per-trial via changeDynamics (pybullet_sweep.py's run_push
uses a module-level constant PUSHER_FRICTION fixed for the whole sweep).
Physics setup (board/pusher dimensions, push speed/distance, classifier) is
identical to pybullet_sweep.py for a like-for-like comparison.
"""
from __future__ import annotations
import numpy as np

BOARD_HALF = (0.103, 0.0171, 0.0125)
PUSH_DIST = 0.040
BOARD_MASS = 1.2
TABLE_FRICTION = 0.30
CONTACT_X = -BOARD_HALF[0]  # push at the board's edge, matching ContactGeometry(p0=(CONTACT_X, 0))


def build_world(gui=False):
    import pybullet as p
    import pybullet_data
    p.connect(p.GUI if gui else p.DIRECT)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.setGravity(0, 0, -9.81)
    p.loadURDF("plane.urdf")
    col = p.createCollisionShape(p.GEOM_BOX, halfExtents=BOARD_HALF)
    vis = p.createVisualShape(p.GEOM_BOX, halfExtents=BOARD_HALF, rgbaColor=[0.76, 0.64, 0.42, 1])
    board = p.createMultiBody(BOARD_MASS, col, vis, basePosition=[0, 0, BOARD_HALF[2] + 1e-3])
    p.changeDynamics(board, -1, lateralFriction=TABLE_FRICTION)
    pcol = p.createCollisionShape(p.GEOM_CYLINDER, radius=0.005, height=0.03)
    pusher = p.createMultiBody(0, pcol, basePosition=[CONTACT_X, 0, 0.015])
    return None, board, pusher


def run_push_mu(board, pusher, angle_deg, mu, steps=2400):
    """Same mechanics/classifier as pybullet_sweep.run_push, but mu is set
    per-call via changeDynamics so a single world can be reused across a
    box-spanning sweep instead of a single fixed PUSHER_FRICTION."""
    import pybullet as p
    p.changeDynamics(pusher, -1, lateralFriction=float(mu))
    p.resetBasePositionAndOrientation(board, [0, 0, BOARD_HALF[2] + 1e-3], [0, 0, 0, 1])
    p.resetBaseVelocity(board, [0, 0, 0], [0, 0, 0])
    start = np.array([CONTACT_X, 0.0, 0.015])
    p.resetBasePositionAndOrientation(pusher, start, [0, 0, 0, 1])
    a = np.deg2rad(angle_deg)
    direction = np.array([np.cos(a), np.sin(a), 0.0])
    dt = 1.0 / 240
    p.setTimeStep(dt)
    pos0, _ = p.getBasePositionAndOrientation(board)
    per_step = PUSH_DIST / steps
    push_speed = per_step / dt
    p0_body = np.array([CONTACT_X, 0.0, 0.0])
    n_ss = max(1, steps // 5)
    vrel_t_samples = []
    for i in range(1, steps + 1):
        new = start + direction * (per_step * i)
        p.resetBasePositionAndOrientation(pusher, new.tolist(), [0, 0, 0, 1])
        p.resetBaseVelocity(pusher, [direction[0] * push_speed, direction[1] * push_speed, 0.0], [0, 0, 0])
        p.stepSimulation()
        if i > steps - n_ss:
            vlin, vang = p.getBaseVelocity(board)
            _, orn = p.getBasePositionAndOrientation(board)
            yaw = p.getEulerFromQuaternion(orn)[2]
            v_com = np.array(vlin[:2])
            wz = vang[2]
            cy, sy = np.cos(yaw), np.sin(yaw)
            r_world = np.array([cy * p0_body[0] - sy * p0_body[1], sy * p0_body[0] + cy * p0_body[1]])
            v_cp = v_com + np.array([-wz * r_world[1], wz * r_world[0]])
            t_world = np.array([-sy, cy])
            pusher_v = np.array([direction[0] * push_speed, direction[1] * push_speed])
            vrel_t_samples.append(float(np.dot(pusher_v - v_cp, t_world)))
    pos1, _ = p.getBasePositionAndOrientation(board)
    disp = np.array(pos1[:2]) - np.array(pos0[:2])
    along = float(np.dot(disp, direction[:2]))
    eff = 100 * max(0.0, along) / PUSH_DIST
    vrel_t_ss = float(np.mean(np.abs(vrel_t_samples))) if vrel_t_samples else 0.0
    tol = 0.15 * push_speed
    mode = "stick" if vrel_t_ss < tol else "slide"
    return mode, eff
