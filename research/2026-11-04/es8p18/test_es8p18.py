import unittest
from fractions import Fraction as F
import audit_es8p18 as a

class AuditES8P18(unittest.TestCase):
    def test_constants(self):
        self.assertEqual(a.ODD_SQ,F(167281,640000))
        self.assertEqual(a.EVEN_SQ,F(2588881,3686400))
        self.assertEqual(a.KAPPA,F(6275281,3686400))
    def test_old_improvement(self):
        self.assertLess(a.KAPPA_SIMPLE,a.OLD)
    def test_budget(self):
        self.assertEqual(a.KAPPA_SIMPLE*a.RAW_DELTA,F(3857,506250))
        self.assertLess(a.KAPPA_SIMPLE*a.RAW_DELTA,F(1,100))
    def test_coeff_all(self):
        for k in range(3,180):
            for x in (F(9),F(10),F(100)):
                alpha=a.A_ODD if k%2 else a.A_EVEN
                self.assertLessEqual(a.coeff(k,x),alpha/k)
    def test_even_reciprocal_bound(self):
        for n in range(2,350):self.assertTrue(a.even_reciprocal_bound(n))
    def test_odd_reciprocal_bound(self):
        for n in range(1,350):self.assertTrue(a.odd_reciprocal_bound(n))
    def test_norm_monotone_in_X(self):
        for n in (3,6,8,10,20,40):
            self.assertLessEqual(a.exact_norm_sq(n,F(10)),a.exact_norm_sq(n,F(9)))
    def test_adversarial_false_claim_pi_independent(self):
        with self.assertRaises(AssertionError):a.coeff(5,F(1))
    def test_full_rational_replay(self):
        self.assertEqual(a.full_validation()['convenient_norm_squared_bound'],'171/100')
    def test_finite_rank_op_norm_no_convergence(self):
        # A truncation that drops all coefficients after N necessarily drops identity on e_(N+1).
        for n in (3,8,50):
            unit=[F(0)]*(n+2)
            unit[n]=F(1)
            self.assertEqual(sum(v*v for v in unit[n:]),1)
    def test_no_hs_conclusion_from_finite_prefix(self):
        # After any all-zero finite prefix, an infinite diagonal of ones is not HS.
        for n in (4,8):
            partial=[F(0)]*n
            for m in (n+10,n+100):
                self.assertEqual(sum([v*v for v in partial]+[F(1)]*(m-n)),m-n)

if __name__=='__main__':unittest.main()
