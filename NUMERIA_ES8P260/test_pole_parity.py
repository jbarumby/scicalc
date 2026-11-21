"""NUMERIA P260 pole normalization identities in exact fractions."""
import unittest
from fractions import Fraction as F
from .pole_parity import legendre_coefficients,legendre_moment,modified_spherical_bessel_coefficient,audit

class PoleNormalizationTests(unittest.TestCase):
    def test_standard_legendre_polynomials(self):
        self.assertEqual(legendre_coefficients(0),(F(1),))
        self.assertEqual(legendre_coefficients(1),(F(0),F(1)))
        self.assertEqual(legendre_coefficients(2),(F(-1,2),F(0),F(3,2)))
    def test_zero_moments_below_degree(self):
        for n in range(1,11):
            for k in range(n):
                self.assertEqual(legendre_moment(n,k),F(0))
    def test_parity_moments(self):
        for n in range(10):
            for k in range(2*n+4):
                if (k-n)%2:
                    self.assertEqual(legendre_moment(n,k),F(0))
    def test_first_bessel_coefficients(self):
        self.assertEqual(modified_spherical_bessel_coefficient(0,0),F(1))
        self.assertEqual(modified_spherical_bessel_coefficient(0,1),F(1,6))
        self.assertEqual(modified_spherical_bessel_coefficient(1,0),F(1,3))
    def test_240_exact_checks(self):
        self.assertEqual(audit()['exact_series_coefficients_compared'],240)
    def test_fail_closed_invalid_parameters(self):
        with self.assertRaises(ValueError):audit(0,10)
        with self.assertRaises(ValueError):legendre_coefficients(-1)
        with self.assertRaises(ValueError):modified_spherical_bessel_coefficient(-1,2)
    def test_scope_flags(self):
        x=audit(5,3)
        self.assertFalse(x['entire_Weil_explicit_formula_pole_identification_proved_here'])
        self.assertFalse(x['Riemann_hypothesis_proved'])
        self.assertEqual((x['corresponding_even_matrix_factor'],x['corresponding_odd_matrix_factor']),('+2','-2'))

if __name__=='__main__':unittest.main()
