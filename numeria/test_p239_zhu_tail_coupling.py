import json
import unittest
from pathlib import Path

from numeria.p237.p239_zhu_tail_coupling import build


ROOT = Path(__file__).resolve().parents[0]


class TestP239ZhuTailCoupling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = build(ROOT / 'p237' / 'data', ROOT / 'P238_SCHUR_MARGIN.json', 512)

    def test_both_parities_have_certified_tails_and_couplings(self):
        self.assertTrue(self.report['finite_to_infinite_tail_and_coupling_certified'])
        sectors = {s['parity']: s for s in self.report['sector_bounds']}
        self.assertEqual(set(sectors), {'even', 'odd'})
        for result in sectors.values():
            self.assertTrue(result['ratio_below_3_over_40'])
            self.assertTrue(result['epsilon_B_below_1e-90'])
            self.assertTrue(result['epsilon_D_below_1e-190'])

    def test_global_window_bound_is_conservative_and_not_rh(self):
        self.assertEqual(self.report['finite_block_margin_lower_exact'], '1/1000000000000000000')
        self.assertEqual(self.report['full_R_lower_bound_strict'], '1/2000000000000000000')
        self.assertTrue(self.report['Q_ge_R_from_one_stroke_envelope'])
        self.assertTrue(self.report['complex_extension_by_parity_lemma'])
        self.assertFalse(self.report['infinite_Weil_positivity_global'])
        self.assertFalse(self.report['RH_proved'])

    def test_source_hashes_are_bound(self):
        hashes = self.report['P237_input_sha256']
        for key in ('matrix', 'center', 'shadow', 'C_even', 'C_odd'):
            self.assertEqual(len(hashes[key]), 64)


if __name__ == '__main__':
    unittest.main()
