import json
import unittest
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent

class P230FiniteSPDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = json.loads((ROOT / 'P230_FINITE_SPD.json').read_text())

    def test_two_sectors_have_independent_center_spd_certificates(self):
        self.assertEqual(self.report['stage'], 'NUMERIA ES8-P230 finite 200x200 SPD perturbation certificate')
        self.assertEqual(self.report['arb_precision_bits'], 512)
        self.assertEqual([s['sector'] for s in self.report['sectors']], ['even', 'odd'])
        for sector in self.report['sectors']:
            cert = sector['center_spd_certificate']
            self.assertEqual(sector['dimension'], 200)
            self.assertTrue(sector['P193_center_exactly_symmetric'])
            self.assertEqual(cert['ldl_positive_pivots'], 200)
            self.assertEqual(cert['lambda_min_center_lower_exact'], '1/100000000000000000000')
            inv_norm_sq = cert['inverse_frobenius_squared_arb_enclosure'].strip('[]').split(' +/- ')
            self.assertLess(Decimal(inv_norm_sq[0]) + Decimal(inv_norm_sq[1]), Decimal(10)**40)
            self.assertEqual(sector['P227_interval_entries_compared'], 20100)
            self.assertIn('symmetric', sector['P227_target_symmetry_basis'])
            delta = F(sector['P227_to_P193center_infinity_norm_upper_exact'])
            lower = F(sector['lambda_min_P227_finite_lower_exact'])
            simple_lower = F(sector['simple_uniform_lower_bound_exact'])
            self.assertGreater(lower, simple_lower)
            self.assertLess(delta, F(1, 200000000000000000000))
            self.assertTrue(sector['simple_uniform_lower_bound_holds'])

    def test_global_claim_is_explicitly_disclaimed(self):
        scope = self.report['scope'].lower()
        self.assertIn('no p207 tail/coupling proof', scope)
        self.assertIn('global weil positivity', scope)
        self.assertIn('rh proof', scope)

if __name__ == '__main__':
    unittest.main()
