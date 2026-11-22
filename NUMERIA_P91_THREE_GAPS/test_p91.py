"""P91 regressions; proof is the full coefficient check in certificate.certify()."""
import unittest
from fractions import Fraction as F
from .certificate import (
    TRIANGLES, SCHEMA, ENVELOPES, RANGES, certify, in_cell,
    triangle_area, envelope, middle, homogeneous_numerator, evaluate,
)

class ExactThreeGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ledger = certify()

    def test_proof_completeness_and_scope(self):
        d = self.ledger
        self.assertEqual(d["schema"], SCHEMA)
        self.assertTrue(d["three_gap_budget_proved"])
        self.assertTrue(d["reflection_reduces_to_a_le_c"])
        self.assertTrue(d["all_coefficients_nonnegative"])
        self.assertEqual(d["nonzero_coefficients_total"], 447)
        self.assertEqual(len(d["triangles"]), 9)
        for flag in ("four_or_more_gaps_proved", "full_Preissmann_Lemma5_reconstructed",
                     "full_Preissmann_Lemma6_reproved", "Riemann_Hypothesis_proved"):
            self.assertFalse(d[flag], flag)

    def test_independent_area_and_cell_membership(self):
        area = F(0)
        for order,left,right,raw in TRIANGLES:
            points = tuple(tuple(map(F,p)) for p in raw)
            self.assertTrue(all(in_cell(a,c,order,left,right) for a,c in points))
            area += triangle_area(points)
            a = sum(v[0] for v in points)/3
            c = sum(v[1] for v in points)/3
            self.assertTrue(in_cell(a,c,order,left,right))
        self.assertEqual(area, F(1,4))

    def test_rebuild_all_rational_coefficients(self):
        expected_counts = [53,42,41,44,54,55,49,54,55]
        for idx, (order,left,right,raw) in enumerate(TRIANGLES):
            vertices = tuple(tuple(map(F,p)) for p in raw)
            p = homogeneous_numerator(order,left,right,vertices)
            self.assertEqual(len(p), expected_counts[idx])
            self.assertTrue(all(val > 0 for val in p.values()))
            self.assertTrue(all(sum(exp)==9 for exp in p))
            expected = self.ledger["triangles"][idx]["coefficients"]
            actual = {f"{i},{j},{k}":str(val) for (i,j,k),val in sorted(p.items())}
            self.assertEqual(actual, expected)

    def test_independent_homogeneous_identity(self):
        weights = ((1,1,1), (2,3,5), (5,3,2), (1,7,4))
        for order,left,right,raw in TRIANGLES:
            vertices = tuple(tuple(map(F,p)) for p in raw)
            p = homogeneous_numerator(order,left,right,vertices)
            for t in weights:
                H = sum(t)
                a = sum(F(t[i])*vertices[i][0] for i in range(3))/H
                c = sum(F(t[i])*vertices[i][1] for i in range(3))/H
                expected = (4*a*a*c*c*(1-a)**2*(1-c)**2
                            *(envelope(a)+envelope(c)-middle(a,c)))
                self.assertEqual(evaluate(p,t), H**9*expected)
                self.assertGreaterEqual(expected,0)

    def test_uniform_equality_and_chord_joints(self):
        self.assertEqual(envelope(F(1,3)), F(27,4))
        self.assertEqual(middle(F(1,3),F(1,3)), F(27,2))
        self.assertEqual(envelope(F(1,3))*2,middle(F(1,3),F(1,3)))
        self.assertEqual(envelope(F(1,2)),F(4))
        self.assertEqual(envelope(F(1)),F(0))
        for key,interval in RANGES.items():
            a,b = interval
            for s in (a,b):
                if s:
                    alpha,beta = ENVELOPES[key]
                    self.assertEqual(envelope(s),alpha/s+beta-3*s)

    def test_finite_rational_mesh_is_only_regression(self):
        for denominator in range(3,19):
            for i in range(1,denominator-1):
                for j in range(1,denominator-i):
                    a,c = F(i,denominator),F(j,denominator)
                    b = 1-a-c
                    self.assertGreater(b,0)
                    self.assertGreaterEqual(envelope(a)+envelope(c)-middle(a,c),0)

    def test_guards(self):
        for s in (F(0),F(-1),F(2)):
            with self.assertRaises(ValueError): envelope(s)
        for a,c in ((F(0),F(1,3)), (F(1,2),F(1,2)), (F(1),F(1))):
            with self.assertRaises(ValueError): middle(a,c)

if __name__=='__main__':
    unittest.main()
