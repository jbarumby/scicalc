"""P92F exact rational regression suite; universal theorem in README.md."""
import unittest
from fractions import Fraction as F
from .certificate import (
    SCHEMA, ZETA2_RATIONAL_UPPER, PI_RATIONAL_LOWER,
    SQRT_6_OVER_5_RATIONAL_LOWER, square_root_exact,
    points_deltas_from_gaps, weighted_schur_rows, delta_matrix,
    check_class_shells, exact_positive_definite_ldlt, certify,
)

class P92FTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ledger=certify()

    def test_census_and_scope(self):
        x=self.ledger
        self.assertEqual(x['schema'],SCHEMA)
        self.assertEqual(x['regression_count'],368)
        self.assertEqual(x['ldlt_certified_count'],40)
        self.assertTrue(x['three_scales_imply_P88_residual_bound'])
        self.assertFalse(x['arbitrarily_many_scales_proved'])
        self.assertFalse(x['full_Preissmann_Lemma5_or_6_proved'])
        self.assertFalse(x['Riemann_Hypothesis_proved'])

    def test_gap_nearest_neighbor_rule(self):
        x,d=points_deltas_from_gaps((1,4,4,9,9))
        self.assertEqual(x,[F(t) for t in (0,1,5,9,18,27)])
        self.assertEqual(d,[F(t) for t in (1,1,4,4,9,9)])
        self.assertEqual(len(set(d)),3)

    def test_delta_vanishes_on_same_scale(self):
        x,d=points_deltas_from_gaps((1,4,4,9,9))
        r=[square_root_exact(v) for v in d]
        mat=delta_matrix(x,d,r)
        for i in range(len(d)):
            for j in range(len(d)):
                if d[i]==d[j]: self.assertEqual(mat[i][j],0)
                if i!=j:self.assertEqual(mat[i][j],mat[j][i])

    def test_two_scale_and_three_scale_fraction_shells(self):
        for gaps in ((1,1,4,4),(1,4,4,9,9),(9,9,4,4,1),(1,4,9,4,9,1)):
            x,d=points_deltas_from_gaps(gaps)
            roots=[square_root_exact(v) for v in d]
            check_class_shells(x,d,roots)
            rows=weighted_schur_rows(x,d,roots)
            self.assertLessEqual(max(rows),(len(set(d))-1)*ZETA2_RATIONAL_UPPER)

    def test_boundary_constant_rational_inequalities(self):
        self.assertLess(F(109,100)**2,F(6,5))
        self.assertLess(F(25,8),F(22,7))
        self.assertLess(F(22,7)**2,10)
        self.assertLess(2*F(5,3),PI_RATIONAL_LOWER**2/3*SQRT_6_OVER_5_RATIONAL_LOWER)

    def test_elementary_scale_reduction_polynomial(self):
        for den in (1,2,3,10,37):
            for i in range(den+1):
                u=F(i,den)
                self.assertLessEqual((1-u)**2*(1+u*u),1)
                self.assertGreater(1-u+2*u*u,0)

    def test_exact_ldlt_independent(self):
        for gaps in ((1,4,4,9,9),(9,4,1,4,9),(1,4,9,9,4,1)):
            x,d=points_deltas_from_gaps(gaps)
            r=[square_root_exact(v) for v in d]
            mat=delta_matrix(x,d,r)
            pivot=exact_positive_definite_ldlt((len(set(d))-1)*F(5,3)+F(1,100),mat)
            self.assertGreater(pivot,0)

    def test_invalid_inputs(self):
        for bad in ((),(1,0,4),(-1,4)):
            with self.assertRaises(ValueError): points_deltas_from_gaps(bad)
        with self.assertRaises(ValueError):square_root_exact(F(2))
        with self.assertRaises(ValueError):square_root_exact(F(-1))

if __name__=='__main__':
    unittest.main()
