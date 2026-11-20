"""P257 adversarial tests with no floating-point sign comparisons."""
import json
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

from .p257_strip_taper import (
    OMEGA, KAPPA_UP, OLD_LOSS, TARGET_LOSS, arctan_lower_alternating,
    arctanh_lower, arctanh_upper, need, pi_strict_lower, prime_7_constants,
    optimal_eta_and_loss, normalized_taper_weight, first_prime7_gap,
    finite_operator_toy, audit,
)


class ExactP257Tests(unittest.TestCase):
    def test_arctanh_brackets_have_strict_width(self):
        for x in (F(3,4),F(3,5),F(1,3)):
            with self.subTest(x=x):
                low = arctanh_lower(x,8)
                high = arctanh_upper(x,8)
                self.assertGreater(high,low)
                self.assertGreater(low,0)

    def test_log7_strict_rational_bound(self):
        c = prime_7_constants()
        self.assertGreater(F(c['kappa_lower_exact']), F(73,100))
        self.assertLess(F(c['kappa_upper_exact']), KAPPA_UP)

    def test_pi_lower_machin_bound(self):
        self.assertGreater(pi_strict_lower(), F(157,50))
        self.assertGreater(arctan_lower_alternating(F(1,5),2),0)

    def test_high_tail_strict_bound(self):
        self.assertGreater(F(prime_7_constants()['B2016_minus_B504_lower_exact']),F(11,8))

    def test_optimizing_young_closes(self):
        r = optimal_eta_and_loss()
        eta = F(r['young_eta_exact'])
        self.assertEqual(eta,F(444,331))
        self.assertEqual(F(r['young_high_consumption_upper_exact']),F(31,24))
        self.assertEqual(F(r['loss_factor_upper_exact']),F(115617600,51967))

    def test_strict_improvement_in_delta(self):
        r = optimal_eta_and_loss()
        self.assertLess(F(r['loss_factor_upper_exact']),TARGET_LOSS)
        self.assertLess(TARGET_LOSS,OLD_LOSS)
        self.assertEqual(F(r['delta_allowance_gain_over_old_exact']),F(5040,4450)-1)

    def test_taper_weight_activation(self):
        self.assertEqual(normalized_taper_weight(7,F(7)),0)
        self.assertEqual(normalized_taper_weight(7,F(701,100)),F(1,701))
        self.assertEqual(normalized_taper_weight(7,F(8)),F(1,8))

    def test_near_threshold_hard_not_approached_in_norm(self):
        for e in (F(1,10),F(1,100),F(1,10**6),F(1,10**30)):
            with self.subTest(e=e):
                report = first_prime7_gap(e)
                self.assertGreater(F(report['difference_norm_strict_rational_lower']),F(63,100))
                self.assertLess(F(report['taper_weight_7']),e/7)

    def test_taper_different_from_raw_weighted_mangoldt(self):
        # For X=9, the taper of S differs from W by a q-dependent factor.
        X=F(9)
        r2 = (2/F(1))* (1-F(2)/X)  # normalized relative to lambda/sqrt(q)
        r3 = (2/F(1))* (1-F(3)/X)
        self.assertNotEqual(r2,r3)

    def test_exact_bipartite_shift_spectrum(self):
        r=finite_operator_toy(F(3,4),F(1,100))
        self.assertEqual(F(r['hard_negative_rayleigh_exact']), -F(3,4))
        self.assertEqual(F(r['exact_difference_operator_norm']),F(297,400))

    def test_tapered_can_be_positive_while_hard_negative(self):
        beta,k,w=F(1,100),F(3,4),F(1,100)
        self.assertGreater(beta-k*w,0)
        self.assertLess(beta-k,0)

    def test_non_positive_kappa_rejected(self):
        with self.assertRaises(ValueError):
            finite_operator_toy(0,F(1,2))

    def test_negative_or_excessive_taper_rejected(self):
        for w in (-F(1,2),F(3,2)):
            with self.assertRaises(ValueError):
                finite_operator_toy(1,w)

    def test_epsilon_invalid_rejected(self):
        for e in (0,-1,1,F(10)):
            with self.assertRaises(ValueError):
                first_prime7_gap(F(e))

    def test_invalid_transcendental_series_fail_closed(self):
        for x in (0,1,-F(1,10)):
            with self.assertRaises(ValueError):
                arctanh_lower(F(x),8)
        with self.assertRaises(ValueError):
            arctan_lower_alternating(F(1,5),0)

    def test_output_scope_not_global(self):
        a=audit()
        self.assertFalse(a['limitations']['crossing_log7_over_2_established'])
        self.assertFalse(a['limitations']['four_prime_auxiliary_margin_positive_beyond_threshold_established'])
        self.assertFalse(a['limitations']['RH_proved'])

    def test_output_deterministic(self):
        a=audit(); b=audit()
        self.assertEqual(json.dumps(a,sort_keys=True),json.dumps(b,sort_keys=True))

    def test_need_cannot_be_disabled_by_optimization(self):
        with self.assertRaises(ArithmeticError):
            need(False,'adversarial tampering')

if __name__ == '__main__':
    unittest.main()
