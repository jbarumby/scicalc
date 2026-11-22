"""P92D exact rational regression suite. Universal proof is in README.md."""
import unittest
from fractions import Fraction as Q
from .certificate import (
    SCHEMA, certify, star_norm_square, square_sum_upper,
    bound_via_integral, star_weights, exact_matrix_sample
)

class P92DTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate=certify()

    def test_scope_and_56_exact_checks(self):
        d=self.certificate
        self.assertEqual(d['schema'],SCHEMA)
        self.assertEqual(d['rational_regression_count'],56)
        self.assertTrue(d['star_family_proved_analytically'])
        for key in ('global_Delta_bound_proved','full_Preissmann_Lemma6_proved',
                    'Riemann_Hypothesis_proved'):
            self.assertFalse(d[key],key)

    def test_uniform_gap_delta_zero(self):
        for m in (2,3,10,100):
            self.assertEqual(star_norm_square(1,m),0)
            self.assertTrue(all(w==0 for w in star_weights(1,m)))

    def test_rank_two_spectral_formula(self):
        for s in (Q(3,2),Q(2),Q(3),Q(10)):
            for m in (2,3,7):
                self.assertEqual(exact_matrix_sample(s,m),star_norm_square(s,m))

    def test_decreasing_series_integral_ceiling(self):
        for s in (Q(1),Q(11,10),Q(3,2),Q(2),Q(5),Q(13)):
            la=s*s
            for m in (2,10,100):
                finite=sum((1/(la+k)**4 for k in range(m)),Q(0))
                self.assertLessEqual(finite,square_sum_upper(s))
                self.assertLessEqual(star_norm_square(s,m),bound_via_integral(s))

    def test_elementary_polynomial_decreasing(self):
        for denominator in (1,2,3,7,19,100):
            for n in range(denominator+1):
                u=Q(n,denominator)
                self.assertLessEqual((1-u)**4*(1+3*u*u),1)
                self.assertGreater(4-6*u+18*u*u,0)

    def test_strictly_below_one_twelfth_for_finite_scale(self):
        for s in (Q(3,2),Q(2),Q(3),Q(50)):
            for m in (2,100):
                self.assertLess(star_norm_square(s,m),Q(1,12))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): star_norm_square(Q(1,2),3)
        with self.assertRaises(ValueError): star_weights(1,1)
        with self.assertRaises(ValueError): square_sum_upper(Q(1,3))


if __name__=='__main__':
    unittest.main()
