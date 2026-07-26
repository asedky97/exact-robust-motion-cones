"""
Exact Robust Motion Cones for planar quasi-static pushing under
parametric and geometric contact uncertainty.

Model (ellipsoidal limit surface, Lynch & Mason 1996; Goyal et al. 1991):
  - Object twist direction for contact force f applied at body-frame point p:
        v = (f_x, f_y, s * tau),   s = 1/c^2,  tau = p x f = p_x f_y - p_y f_x
    where c = Tbar/Fbar is the characteristic length of the limit surface.
  - Contact-point velocity:  u = f + s * (perp(p) . f) * perp(p) = M f,
        M = I + s * perp(p) perp(p)^T   (symmetric positive definite)
  - Friction cone edges at the pusher contact (friction mu, inward normal n):
        f_pm = R(+/- gamma) n,  gamma = atan(mu)
  - Motion cone MC(theta) = cone{ u_-, u_+ } in pusher-velocity space:
    commanding a pusher velocity whose direction lies inside MC yields sticking.

Uncertain parameters theta = (mu, c, alpha, delta):
  mu    : pusher-object friction coefficient
  c     : limit-surface characteristic length (support pressure distribution)
  alpha : rotation of the contact normal (surface roughness / normal uncertainty)
  delta : contact position offset along the nominal edge tangent (geometric)

Key closed forms used by the exact computation (Theorems 1-2 of the paper):
  T1 monotonicity:
    d(phi)/d(mu)  : edges rotate outward => lower edge max and upper edge min at mu_lo
    d(phi)/d(alpha) = det(M) |f|^2 / |M f|^2  > 0  => both edges increase with alpha
    d(phi)/d(s)   = tau(f) (f . p) / |u|^2        => sign constant in s (vertex in c)
  T2 position uncertainty:
    u(delta) = A + B delta + C delta^2  (componentwise quadratic), and
    u x u'   = (A x B) + 2 (A x C) delta + (B x C) delta^2
    => at most two interior critical points, roots of an explicit quadratic.

Author: Ahmed S. M. Elsayed and Sotirios Grammatikos -- code drafted with AI assistance, verify before use.
"""

from __future__ import annotations
import itertools
import numpy as np
from dataclasses import dataclass, field

# ----------------------------------------------------------------------------- helpers

def perp(a: np.ndarray) -> np.ndarray:
    """90-degree CCW rotation: perp((x, y)) = (-y, x)."""
    return np.array([-a[1], a[0]])

def cross2(a: np.ndarray, b: np.ndarray) -> float:
    """2D scalar cross product a x b."""
    return a[0] * b[1] - a[1] * b[0]

def rot(angle: float) -> np.ndarray:
    ca, sa = np.cos(angle), np.sin(angle)
    return np.array([[ca, -sa], [sa, ca]])

def ang(v: np.ndarray) -> float:
    return float(np.arctan2(v[1], v[0]))

def wrap_to(angle: float, ref: float) -> float:
    """Wrap angle into (ref - pi, ref + pi]."""
    return ref + np.remainder(angle - ref + np.pi, 2 * np.pi) - np.pi


# ----------------------------------------------------------------------------- model

@dataclass
class ContactGeometry:
    """Nominal contact on the object boundary, body frame, CoM at origin."""
    p0: np.ndarray          # nominal contact point
    n0: np.ndarray          # nominal inward unit normal
    t0: np.ndarray = None   # edge tangent (unit); default = perp(n0)

    def __post_init__(self):
        self.p0 = np.asarray(self.p0, dtype=float)
        self.n0 = np.asarray(self.n0, dtype=float)
        self.n0 = self.n0 / np.linalg.norm(self.n0)
        if self.t0 is None:
            self.t0 = perp(self.n0)
        else:
            self.t0 = np.asarray(self.t0, dtype=float)
            self.t0 = self.t0 / np.linalg.norm(self.t0)

    def contact(self, alpha: float, delta: float):
        """Realized contact point and normal for geometric parameters (alpha, delta)."""
        p = self.p0 + delta * self.t0
        n = rot(alpha) @ self.n0
        return p, n


@dataclass
class UncertaintyBox:
    mu: tuple[float, float]
    c: tuple[float, float]
    alpha: tuple[float, float] = (0.0, 0.0)
    delta: tuple[float, float] = (0.0, 0.0)

    def vertices(self):
        return list(itertools.product(self.mu, self.c, self.alpha, self.delta))

    def nominal(self):
        return tuple(0.5 * (lo + hi) for (lo, hi) in
                     (self.mu, self.c, self.alpha, self.delta))


def friction_edges(n: np.ndarray, mu: float):
    """Friction cone edge directions f_-, f_+ (unit), gamma = atan(mu)."""
    g = np.arctan(mu)
    return rot(-g) @ n, rot(+g) @ n

def push_map(p: np.ndarray, c: float) -> np.ndarray:
    """M = I + s perp(p) perp(p)^T mapping contact force direction to contact velocity."""
    s = 1.0 / c**2
    q = perp(p)
    return np.eye(2) + s * np.outer(q, q)

def cone_rays(geom: ContactGeometry, mu: float, c: float,
              alpha: float, delta: float):
    """Boundary rays (u_minus, u_plus) of the sticking motion cone, plus interior ref."""
    p, n = geom.contact(alpha, delta)
    fm, fp = friction_edges(n, mu)
    M = push_map(p, c)
    return M @ fm, M @ fp, M @ n   # u_-, u_+, interior reference direction

def cone_angles(geom, mu, c, alpha, delta, ref: float):
    """Edge angles (phi_minus, phi_plus) wrapped near a common reference angle."""
    um, up, ui = cone_rays(geom, mu, c, alpha, delta)
    a_m, a_p = wrap_to(ang(um), ref), wrap_to(ang(up), ref)
    lo, hi = min(a_m, a_p), max(a_m, a_p)
    return lo, hi


# ------------------------------------------------------- robust cone: brute force

def robust_cone_bruteforce(geom: ContactGeometry, box: UncertaintyBox,
                           n_grid: int = 15):
    """Sampled approximation of the infinite intersection (baseline, Liang et al. style
    but on a dense grid rather than vertices only). Returns (Phi_lo, Phi_hi, ref)."""
    mu0, c0, a0, d0 = box.nominal()
    _, _, ui = cone_rays(geom, mu0, c0, a0, d0)
    ref = ang(ui)
    lo_best, hi_best = -np.inf, np.inf
    grids = [np.linspace(lo, hi, n_grid) if hi > lo else np.array([lo])
             for (lo, hi) in (box.mu, box.c, box.alpha, box.delta)]
    for mu in grids[0]:
        for c in grids[1]:
            for al in grids[2]:
                for de in grids[3]:
                    lo, hi = cone_angles(geom, mu, c, al, de, ref)
                    lo_best, hi_best = max(lo_best, lo), min(hi_best, hi)
    return lo_best, hi_best, ref


# ------------------------------------------------------- robust cone: exact (T1+T2)

def _delta_critical_points(geom: ContactGeometry, mu: float, c: float,
                           alpha: float, which: str):
    """Interior critical points of phi_edge(delta): roots of
       (A x B) + 2 (A x C) delta + (B x C) delta^2 = 0,
    with u(delta) = A + B delta + C delta^2 (Theorem 2)."""
    s = 1.0 / c**2
    n = rot(alpha) @ geom.n0
    fm, fp = friction_edges(n, mu)
    f = fm if which == "minus" else fp
    q1, q2 = perp(geom.p0), perp(geom.t0)
    k0, k1 = float(q1 @ f), float(q2 @ f)
    A = f + s * k0 * q1
    B = s * (k1 * q1 + k0 * q2)
    C = s * k1 * q2
    coeffs = [cross2(B, C), 2.0 * cross2(A, C), cross2(A, B)]  # quadratic in delta
    if abs(coeffs[0]) < 1e-14 and abs(coeffs[1]) < 1e-14:
        return []
    roots = np.roots(coeffs) if abs(coeffs[0]) > 1e-14 else \
        np.array([-coeffs[2] / coeffs[1]])
    return [float(r) for r in roots if abs(r.imag if np.iscomplexobj(r) else 0) < 1e-10]

def robust_cone_exact(geom: ContactGeometry, box: UncertaintyBox):
    """Exact robust motion cone via the finite candidate set of Theorems 1-2.

    Candidate parameter values per edge:
      mu    : mu_lo only (cone edges rotate outward with mu)        [T1-i]
      c     : both endpoints (sign of d(phi)/ds constant in s)      [T1-ii]
      alpha : per-edge: alpha_max for phi-, alpha_min for phi+      [T1-iii]
      delta : endpoints + interior roots of the explicit quadratic  [T2]
    Returns (Phi_lo, Phi_hi, ref, n_evals).
    """
    mu0, c0, a0, d0 = box.nominal()
    _, _, ui = cone_rays(geom, mu0, c0, a0, d0)
    ref = ang(ui)
    mu_lo = box.mu[0]

    def _eval_candidates(ref_angle):
        """Inner loop: evaluate all candidates unwrapped about ref_angle.
        Returns (lo_best, hi_best, n_evals, any_wraparound)."""
        lo, hi, ne = -np.inf, np.inf, 0
        wrap = False
        # Theorem 1(iii): phi- needs sup -> alpha_max, phi+ needs inf -> alpha_min
        for which, sel, alpha_val in (
            ("minus", max, box.alpha[1]), ("plus", min, box.alpha[0]),
        ):
            for c in set(box.c):
                al = alpha_val
                deltas = set(box.delta)
                for r in _delta_critical_points(geom, mu_lo, c, al, which):
                    if box.delta[0] < r < box.delta[1]:
                        deltas.add(r)
                for de in deltas:
                    um, up, _ = cone_rays(geom, mu_lo, c, al, de)
                    # Check if this cone's interval wraps around ref_angle
                    a_m = wrap_to(ang(um), ref_angle)
                    a_p = wrap_to(ang(up), ref_angle)
                    lo_c, hi_c = min(a_m, a_p), max(a_m, a_p)
                    if hi_c - lo_c > np.pi:
                        wrap = True  # ref is not interior to this cone
                    u = um if which == "minus" else up
                    a = wrap_to(ang(u), ref_angle)
                    ne += 1
                    if which == "minus":
                        lo = max(lo, a)
                    else:
                        hi = min(hi, a)
        return lo, hi, ne, wrap

    lo_best, hi_best, n_evals, wrapped = _eval_candidates(ref)

    # If the intersection appears non-empty but some cones did not contain the
    # reference, the unwrapped representation may be spurious (large interval > pi
    # from a cone that does not actually contain the reference). Re-centre once.
    if lo_best < hi_best and wrapped:
        ref = 0.5 * (lo_best + hi_best)
        lo_best, hi_best, n_evals2, wrapped2 = _eval_candidates(ref)
        n_evals += n_evals2

    # If empty (or still empty after re-centering), one more re-centre attempt
    # from a different branch: try the midpoint of the two edge candidates.
    if lo_best >= hi_best:
        mid = 0.5 * (lo_best + hi_best)
        if abs(wrap_to(mid, ref) - mid) > 1e-10:
            # Midpoint is in a different branch from current ref; try it
            ref = mid
            lo_best, hi_best, n_evals2, _ = _eval_candidates(ref)
            n_evals += n_evals2

    return lo_best, hi_best, ref, n_evals


def sticking_omega(p: np.ndarray, c: float, u: np.ndarray) -> float:
    """Closed-form rotation rate under sticking (Theorem 3):
         omega = (p x u) / (c^2 + |p|^2)."""
    return cross2(p, u) / (c**2 + float(p @ p))
