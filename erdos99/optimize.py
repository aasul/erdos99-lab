"""Epigraph SLSQP with analytic derivatives and optional triangle-free branches."""
import time
import numpy as np
from scipy.optimize import minimize
from .geometry import canonical, normalize, pairs, squared_distances, triple_pair_indices, triangle_margin


class DiameterProblem:
    def __init__(self, n, lower=None):
        self.n = n
        self.i, self.j = pairs(n)
        self.m = len(self.i)
        # p0=(0,0), p1=(x,0). Free coordinates are full flattened indices 2,4,5,...
        self.free = np.r_[2, np.arange(4, 2*n)]
        self.lower = np.ones(self.m) if lower is None else np.asarray(lower)

    def unpack(self, z):
        p = np.zeros(2*self.n)
        p[self.free] = z[:-1]
        return p.reshape(-1, 2)

    def initial(self, points):
        p = canonical(points)
        # Inflate to make the branch start feasible, including all witness edges.
        p *= max(1., np.sqrt(np.max(self.lower / squared_distances(p))))
        return np.r_[p.ravel()[self.free], squared_distances(p).max()*(1+1e-10)]

    def constraints(self, z):
        q = squared_distances(self.unpack(z))
        return np.r_[q-self.lower, z[-1]-q]

    def jacobian(self, z):
        p = self.unpack(z)
        delta = 2*(p[self.i]-p[self.j])
        jac = np.zeros((self.m, 2*self.n))
        rows = np.arange(self.m)
        jac[rows, 2*self.i] = delta[:, 0]
        jac[rows, 2*self.i+1] = delta[:, 1]
        jac[rows, 2*self.j] = -delta[:, 0]
        jac[rows, 2*self.j+1] = -delta[:, 1]
        jac = jac[:, self.free]
        return np.block([[jac, np.zeros((self.m, 1))], [-jac, np.ones((self.m, 1))]])


def optimize(points, maxiter=600, margin=0.0, ftol=1e-10):
    start = time.monotonic()
    n = len(points)
    points = normalize(points)
    if margin > 0 and triangle_margin(points) < margin:
        rng = np.random.default_rng(99000+n)
        original = points.copy()
        for k in range(20):
            candidate = normalize(original + rng.normal(0, margin*(5+k), original.shape))
            if triangle_margin(candidate) >= margin:
                points = candidate
                break
        else:
            # Guaranteed feasible branch start, even for an unusually large margin.
            if margin > 1:
                raise ValueError("This search supports margins up to 1")
            points = np.column_stack([np.arange(n), np.zeros(n)])
    lower = np.ones(n*(n-1)//2)
    if margin > 0:
        # Each triple chooses its currently longest edge as a witness.
        # Fixing these choices explores one branch of the disjunctive feasible set.
        tri = triple_pair_indices(n)
        q = squared_distances(points)
        witnesses = tri[np.arange(len(tri)), q[tri].argmax(axis=1)]
        lower[np.unique(witnesses)] = (1+margin)**2
    problem = DiameterProblem(n, lower)
    z0 = problem.initial(points)
    objective_grad = np.r_[np.zeros(len(z0)-1), 1.]
    res = minimize(lambda z: z[-1], z0, jac=lambda z: objective_grad,
                   method="SLSQP", constraints={"type":"ineq", "fun":problem.constraints,
                                                "jac":problem.jacobian},
                   bounds=[(1e-8, None)]+[(None, None)]*(len(z0)-2)+[(1., None)],
                   options={"ftol":ftol, "maxiter":maxiter, "disp":False})
    raw = problem.unpack(res.x)
    raw_residual = float(problem.constraints(res.x).min())
    # Always retain the feasible starting configuration as a fallback.
    candidates = [normalize(problem.unpack(z0))]
    try:
        candidates.append(normalize(raw))
    except ValueError:
        pass
    valid = [p for p in candidates if triangle_margin(p) >= margin-1e-8]
    best = min(valid, key=lambda p: squared_distances(p).max())
    return best, dict(success=bool(res.success), status=int(res.status), message=str(res.message),
                      iterations=int(res.nit), raw_min_residual=raw_residual,
                      raw_min_distance=float(np.sqrt(squared_distances(raw).min())),
                      fallback=bool(best is candidates[0]), elapsed_seconds=time.monotonic()-start)
