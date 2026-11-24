import math
import unittest
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from erdos99.geometry import diagnostics, triangle_margin
from erdos99.starts import known_small, generate
from erdos99.optimize import DiameterProblem, optimize
from erdos99.verify import verify_points


class ResearchTests(unittest.TestCase):
    def test_independent_known_geometry(self):
        values={4:math.sqrt(2),5:(1+math.sqrt(5))/2,6:2*math.sin(2*math.pi/5),
                7:2,8:1/(2*math.sin(math.pi/14))}
        for n,expected in values.items():
            p=known_small(n)
            v=verify_points(p.tolist())
            self.assertTrue(v["exact_saved_feasible"])
            self.assertAlmostEqual(v["normalized_diameter"],expected,places=11)
            self.assertEqual(v["near_unit_triangles"]["1e-6"] > 0,n==7)

    def test_analytic_jacobian(self):
        rng=np.random.default_rng(17)
        p=rng.normal(size=(6,2))
        problem=DiameterProblem(6)
        z=problem.initial(p)
        np.testing.assert_allclose(problem.jacobian(z),approx_derivative(problem.constraints,z),
                                   rtol=1e-6,atol=1e-7)

    def test_scaling_invariance_and_triangle_trap(self):
        p=known_small(7)
        self.assertAlmostEqual(triangle_margin(p),triangle_margin(37*p),places=12)
        inflated=(p*1.00001).tolist()
        v=verify_points(inflated)
        self.assertEqual(v["exact_saved_unit_triangles"],0)
        self.assertEqual(v["near_unit_triangles"]["1e-6"],6)

    def test_verifier_rejects_bad_geometry(self):
        self.assertFalse(verify_points([[0,0],[.9,0],[0,1]])["exact_saved_feasible"])
        with self.assertRaises(ValueError):
            verify_points([[0,0],[0,0]])

    def test_seed_reproduction(self):
        a=list(generate(10,[19])); b=list(generate(10,[19]))
        for x,y in zip(a,b):
            np.testing.assert_array_equal(x[2],y[2])

    def test_optimizer_and_margin(self):
        p=known_small(7)
        q,info=optimize(p,maxiter=200,margin=.01)
        v=verify_points(q.tolist(),.01)
        self.assertTrue(v["exact_saved_feasible"])
        self.assertTrue(v["margin_within_acceptance"])
        q,_=optimize(known_small(4),maxiter=100)
        self.assertAlmostEqual(verify_points(q.tolist())["normalized_diameter"],math.sqrt(2),places=7)


if __name__=="__main__":
    unittest.main()
