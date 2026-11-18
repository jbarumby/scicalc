import json
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p228_recompose_audit import beta_interval, dyad, pole_component
from p229_quadrature_bound import total_error_bound

ROOT = Path(__file__).resolve().parent

class P228AuditTests(unittest.TestCase):
    def test_dyadic_parser_handles_sign_and_exponent(self):
        self.assertEqual(dyad(['3', '-2']), F(3, 4))
        self.assertEqual(dyad(['-5', '3']), F(-40, 1))

    def test_beta_interval_is_narrow(self):
        for T in (150, 200):
            lo, hi = beta_interval(T, 120)
            self.assertLess(lo, hi)
            self.assertLess(hi-lo, F(1, 10**100))

    def test_independent_GL_error_bounds(self):
        self.assertLess(total_error_bound(200), F(3, 10**51))
        self.assertLess(total_error_bound(150), F(3, 10**51))
        p229 = json.loads((ROOT/'P229_BOUND.json').read_text())
        self.assertEqual(p229['stage'], 'NUMERIA ES8-P229 independent GL60 error bounds')
        for sector in p229['sectors']:
            new = F(sector['E_GL_exact'])
            old = F(sector['P193_E_GL_exact'])
            self.assertEqual(new, total_error_bound(sector['T']))
            self.assertGreater(old/new, 10000)

    def test_pole_component_bounds_are_positive(self):
        for n in (0, 1, 5, 199, 398, 399):
            lo, hi = pole_component(n)
            self.assertLessEqual(F(0), lo)
            self.assertLess(lo, hi)
            self.assertLess(F(0), hi)
            self.assertLess(hi-lo, F(1, 10**70))

    def test_recomposed_artifact_result(self):
        report = json.loads((ROOT/'P228_AUDIT.json').read_text())
        self.assertEqual(report['stage'], 'NUMERIA ES8-P228 local rational audit of P227')
        self.assertEqual(len(report['sectors']), 2)
        for sector in report['sectors']:
            self.assertEqual(sector['upper_triangle_entries'], 20100)
            self.assertEqual(sector['entries_intersecting_P193'], 20100)
            self.assertEqual(sector['entries_disjoint'], 0)
            self.assertEqual(len(sector['chunks']), 8 if sector['sector']=='even' else 6)
            self.assertLess(float(sector['max_midpoint_difference_decimal']), 1e-46)
            self.assertLess(F(sector['P227_independent_quadrature_error_exact']), F(3,10**51))
            self.assertGreater(F(sector['P193_quadrature_error_exact']), F(sector['P227_independent_quadrature_error_exact']))
            self.assertLess(F(sector['P227_to_P193center_infinity_norm_upper_exact']), F(1,10**48))

if __name__ == '__main__':
    unittest.main()
