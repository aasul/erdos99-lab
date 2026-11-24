"""Geometry used by the optimizer; verification deliberately lives elsewhere."""
from itertools import combinations
import numpy as np

CONTACT_TOL = 1e-6
SAFETY = 2e-12


def pairs(n):
    return np.triu_indices(n, 1)


def squared_distances(points):
    i, j = pairs(len(points))
    d = np.asarray(points)[i] - np.asarray(points)[j]
    return np.einsum("ij,ij->i", d, d)


def normalize(points, safety=SAFETY):
    p = np.array(points, dtype=float, copy=True)
    if p.ndim != 2 or p.shape[1] != 2 or len(p) < 2 or not np.isfinite(p).all():
        raise ValueError("Expected at least two finite points of shape (n,2)")
    p -= p[0]
    q = squared_distances(p)
    if q.min() <= 1e-24:
        raise ValueError("Coincident or numerically indistinguishable points")
    return p * ((1 + safety) / np.sqrt(q.min()))


def canonical(points):
    p = normalize(points, safety=0)
    x, y = p[1]
    r = np.hypot(x, y)
    p = p @ np.array([[x, -y], [y, x]]) / r
    p[0] = 0
    p[1, 1] = 0
    return p


def triple_pair_indices(n):
    i, j = pairs(n)
    lookup = {tuple(e): k for k, e in enumerate(zip(i, j))}
    return np.array([[lookup[a, b], lookup[a, c], lookup[b, c]]
                     for a, b, c in combinations(range(n), 3)], dtype=int)


def triangle_margin(points):
    """min over triples of (longest side / global minimum distance - 1)."""
    q = squared_distances(points)
    return float(np.sqrt(q[triple_pair_indices(len(points))].max(axis=1).min()
                         / q.min()) - 1)


def diagnostics(points, tol=CONTACT_TOL):
    p = normalize(points, safety=0)
    i, j = pairs(len(p))
    d = np.sqrt(squared_distances(p))
    contacts = d <= 1 + tol
    tri = triple_pair_indices(len(p))
    triangles = int(np.all(contacts[tri], axis=1).sum())
    # Bar-framework rank: contacts treated as equalities, not unilateral constraints.
    rows = np.zeros((int(contacts.sum()), 2 * len(p)))
    for k, (a, b) in enumerate(zip(i[contacts], j[contacts])):
        rows[k, 2*a:2*a+2] = p[a] - p[b]
        rows[k, 2*b:2*b+2] = p[b] - p[a]
    rank = int(np.linalg.matrix_rank(rows, tol=1e-7)) if len(rows) else 0
    active = contacts | (d >= d.max()-tol)
    active_rows = np.zeros((int(active.sum()), 2*len(p)))
    for k,(a,b) in enumerate(zip(i[active],j[active])):
        active_rows[k,2*a:2*a+2] = p[a]-p[b]
        active_rows[k,2*b:2*b+2] = p[b]-p[a]
    active_rank = int(np.linalg.matrix_rank(active_rows,tol=1e-7))
    # Compare rotations about centroid, allowing relabeling using assignment.
    from scipy.optimize import linear_sum_assignment
    c = p - p.mean(axis=0)
    rotations = []
    for order in range(2, 13):
        a = 2*np.pi/order
        rot = c @ np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        costs = np.linalg.norm(rot[:, None] - c[None, :], axis=2)
        aidx, bidx = linear_sum_assignment(costs)
        if costs[aidx, bidx].max() < 1e-5:
            rotations.append(order)
    return dict(contacts=int(contacts.sum()), near_unit_triangles=triangles,
                triangle_margin=triangle_margin(p), contact_rank=rank,
                bar_flex_dimension=2*len(p)-3-rank,
                contact_and_diameter_rank=active_rank,
                fixed_active_length_flex_dimension=2*len(p)-3-active_rank,
                rotational_orders=rotations)
