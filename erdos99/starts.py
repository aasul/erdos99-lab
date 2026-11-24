"""Structured starts, independent random starts, and user supplied geometry."""
import numpy as np
from .geometry import normalize


def polygon(k, radius=1.0, phase=0.0):
    a = np.arange(k)*2*np.pi/k + phase
    return radius*np.column_stack([np.cos(a), np.sin(a)])


def lattice(n, triangular=True, offset=(0., 0.)):
    r = int(np.ceil(np.sqrt(n))) + 2
    ij = np.array([(i, j) for i in range(-r, r+1) for j in range(-r, r+1)])
    p = ij.astype(float)
    if triangular:
        p = np.column_stack([ij[:, 0]+ij[:, 1]/2, ij[:, 1]*np.sqrt(3)/2])
    p -= offset
    order = np.argsort(np.sum(p*p, axis=1), kind="stable")
    return normalize(p[order[:n]])


def known_small(n):
    if n in (4, 5):
        return normalize(polygon(n))
    if n == 6:
        return normalize(np.vstack([[0, 0], polygon(5)]))
    if n == 7:
        return normalize(np.vstack([[0, 0], polygon(6)]))
    if n == 8:
        return normalize(np.vstack([[0, 0], polygon(7)]))
    return None


def literature_nine():
    """Rounded coordinates in Daroczy-Kiss (2005); normalize before use."""
    return normalize([[0,0], [1.0912,0], [2.0912,0], [.5456,.8381],
                      [.5456,-.8381], [1.4404,1.2846], [1.4404,-1.2846],
                      [2.3848,.9559], [2.3848,-.9559]])


def generate(n, seeds):
    yield "triangular", 0, lattice(n)
    yield "square", 0, lattice(n, triangular=False)
    yield "polygon_center", 0, normalize(np.vstack([[0, 0], polygon(n-1)]))
    # Golden-angle spiral supplies a custom non-lattice radial pattern.
    k = np.arange(n)
    a = k*np.pi*(3-np.sqrt(5))
    yield "spiral", 0, normalize(np.sqrt(k+.5)[:, None]*np.column_stack([np.cos(a), np.sin(a)]))
    small = known_small(n)
    if small is not None:
        yield "known_small", 0, small
    if n == 9:
        yield "literature_2005", 0, literature_nine()
    for seed in seeds:
        rng = np.random.default_rng(np.random.SeedSequence([n, seed]))
        yield "random_disk", seed, normalize(np.sqrt(rng.uniform(.01, 1, n))[:, None]
              * np.column_stack([np.cos(a := rng.uniform(0, 2*np.pi, n)), np.sin(a)]))
        yield "jittered_triangular", seed, normalize(lattice(n, offset=tuple(rng.uniform(-.5,.5,2)))
                                                     + rng.normal(0, .10, (n,2)))
