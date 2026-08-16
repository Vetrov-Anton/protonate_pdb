"""Geometry: building an atom from internal coordinates (NeRF) plus metrics."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np


def _v(a) -> np.ndarray:
    if hasattr(a, "x"):
        return np.array([a.x, a.y, a.z], dtype=float)
    return np.asarray(a, dtype=float)


def distance(a, b) -> float:
    return float(np.linalg.norm(_v(a) - _v(b)))


def angle(a, b, c) -> float:
    """Angle a-b-c in degrees."""
    v1, v2 = _v(a) - _v(b), _v(c) - _v(b)
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def dihedral(a, b, c, d) -> float:
    """Dihedral angle a-b-c-d in degrees (IUPAC convention)."""
    p0, p1, p2, p3 = _v(a), _v(b), _v(c), _v(d)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def place_atom(a, b, c, bond: float, ang: float, tors: float) -> np.ndarray:
    """Place atom D so that |C-D| = bond, angle(B,C,D) = ang, dihedral(A,B,C,D) = tors.

    Classic NeRF (natural extension reference frame).
    """
    pa, pb, pc = _v(a), _v(b), _v(c)
    ang_r, tors_r = math.radians(ang), math.radians(tors)

    d2 = np.array(
        [
            -bond * math.cos(ang_r),
            bond * math.cos(tors_r) * math.sin(ang_r),
            bond * math.sin(tors_r) * math.sin(ang_r),
        ]
    )
    bc = pc - pb
    bc /= np.linalg.norm(bc)
    n = np.cross(pb - pa, bc)
    norm = np.linalg.norm(n)
    if norm < 1e-8:  # degenerate case: A, B, C are collinear
        n = np.cross(bc, np.array([1.0, 0.0, 0.0]))
        norm = np.linalg.norm(n)
        if norm < 1e-8:
            n = np.cross(bc, np.array([0.0, 1.0, 0.0]))
            norm = np.linalg.norm(n)
    n /= norm
    m = np.array([bc, np.cross(n, bc), n]).T
    return pc + m.dot(d2)


def tetrahedral_hydrogens(root, neighbor, ref, bond: float = 1.09,
                          ang: float = 109.5, start: float = 60.0) -> list:
    """Three methyl hydrogens on root(-neighbor), staggered relative to ref."""
    return [
        place_atom(ref, neighbor, root, bond, ang, start + 120.0 * i)
        for i in range(3)
    ]


def min_distance_to(point: Sequence[float], cloud: np.ndarray) -> float:
    if len(cloud) == 0:
        return float("inf")
    return float(np.min(np.linalg.norm(cloud - np.asarray(point), axis=1)))


def best_torsion(a, b, c, bond: float, ang: float, cloud: np.ndarray,
                 preferred: float = 180.0, step: float = 15.0):
    """Pick the dihedral that puts the new atom least into its surroundings.

    Starts from `preferred` and returns it right away if it is clear enough,
    otherwise returns the globally best option.
    """
    best, best_score = preferred, -1.0
    for i in range(int(360 / step)):
        tors = preferred + step * i
        pos = place_atom(a, b, c, bond, ang, tors)
        score = min_distance_to(pos, cloud)
        if score > best_score:
            best, best_score = tors, score
        if i == 0 and score > 2.9:      # the original (trans) position is free
            return tors, score
    return best, best_score
