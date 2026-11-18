import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent

class P238SchurBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = json.loads((ROOT / 'P238_SCHUR_MARGIN.json').read_text())

    def test_common_finite_margin_verified_for_both_sectors(self):
        self.assertEqual(self.report['schema'], 'NUMERIA.ES8P238.P237_SCHUR_MARGIN_AUDIT.1')
        self.assertEqual(self.report['common_finite_block_margin_lower_exact'], '1/1000000000000000000')
        for sector in ('even', 'odd'):
            result = self.report['sectors'][sector]
            self.assertEqual(result['dimension'], 377)
            self.assertEqual(result['center_ldl_positive_pivots'], 377)
            self.assertTrue(result['center_spd'])
            self.assertTrue(result['common_margin_strictly_exceeded'])
            self.assertTrue(result['full_interval_box_spd'])
            self.assertEqual(result['common_margin_lower_exact'], '1/1000000000000000000')

    def test_schur_condition_is_conditional_not_global_result(self):
        self.assertEqual(self.report['conditional_normalized_coupling_target_when_epsilon_zero'],
                         'beta/sqrt(eta) < 1/1000000000')
        self.assertFalse(self.report['tail_coercivity_proved'])
        self.assertFalse(self.report['coupling_bound_proved'])
        self.assertFalse(self.report['finite_matrix_identified_with_global_Weil_compression'])
        self.assertFalse(self.report['infinite_Weil_positivity_proved'])
        self.assertFalse(self.report['RH_proved'])

if __name__ == '__main__':
    unittest.main()
