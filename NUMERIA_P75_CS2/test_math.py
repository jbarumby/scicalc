"""Dependency-free exact/rational identity gates for P75 (no flint needed)."""
from fractions import Fraction as F
import math
import unittest
import mpmath as mp

# Original mpmath P75 exploratory script is imported ONLY for independent
# numerical cross-checks; no mpmath output is certified.
from NUMERIA_P75_CS2.preflight import integrated_first_term, coefficient_exponents

class TestFormulaIdentities(unittest.TestCase):
    def test_v2_exponents_at_exact_endpoints(self):
        a,p=coefficient_exponents(F(1,2)); self.assertEqual((a,p),(F(1),F(3)))
        a,p=coefficient_exponents(F(1)); self.assertEqual((a,p),(F(0),F(2)))
        self.assertEqual(coefficient_exponents(F(3,4)),(F(3,5),F(13,5)))

    def test_exact_second_and_first_log_integrals(self):
        # 2 int_H^inf ((log u)^q)/u^3 du, for q=1,2.
        mp.mp.dps=60
        H=mp.mpf(3)*10**12; L=mp.log(H)
        i2=(L*L+L+mp.mpf('0.5'))/H**2
        i1=(L+mp.mpf('0.5'))/H**2
        self.assertLess(abs(i2-2*mp.quad(lambda v:v*v*mp.exp(-2*v),[L,mp.inf])),i2*mp.mpf('1e-45'))
        self.assertLess(abs(i1-2*mp.quad(lambda v:v*mp.exp(-2*v),[L,mp.inf])),i1*mp.mpf('1e-45'))

    def test_integrated_incomplete_gamma_formula(self):
        mp.mp.dps=60
        H=3*10**12
        for sigma in (F(17,32),F(3,4),F(31,32)):
            closed=integrated_first_term(sigma,H)
            a,p=coefficient_exponents(sigma)
            L=mp.log(H)
            alpha=mp.mpf(a.numerator)/a.denominator
            q=mp.mpf(p.numerator)/p.denominator
            direct=2*mp.quad(lambda v:mp.exp(-(2-alpha)*v)*v**q,[L,mp.inf])
            self.assertLess(abs(closed-direct),abs(direct)*mp.mpf('1e-43'))

    def test_secant_quadratic_coefficients(self):
        # x² has f''=2 and exact interpolation error h²*u(1-u).
        for k in range(101):
            u=F(k,100)
            self.assertEqual(u*u, u-u*(1-u))
        # With unequal endpoints, the maximum loss below their minimum
        # is (C-d)²/(4C) when 0 <= d < C.
        for C,d in [(F(1),F(0)),(F(1),F(1,3)),(F(3,2),F(1,2))]:
            u=(C-d)/(2*C)
            loss=C*u*(1-u)-u*d
            self.assertEqual(loss,(C-d)**2/(4*C))

    def test_layer_cake_identity(self):
        mp.mp.dps=50
        # D(s)=sum w_j 1_{theta_j >= s-.5}; exact layer-cake
        for t in [mp.mpf(1),mp.mpf('8.3')]:
            thetas=[mp.mpf('.07'),mp.mpf('.22'),mp.mpf('.42')]
            weights=[mp.mpf('0.4'),mp.mpf('0.2'),mp.mpf('0.1')]
            exact=sum(w*(mp.cosh(th*t)-1) for w,th in zip(weights,thetas))
            integral=mp.quad(lambda x:sum(w for w,th in zip(weights,thetas) if x<th)*t*mp.sinh(x*t),[0,*thetas,mp.mpf('.5')])
            self.assertLess(abs(exact-integral),mp.mpf('1e-45'))

if __name__=='__main__':unittest.main()
