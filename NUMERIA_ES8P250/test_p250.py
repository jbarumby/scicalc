"""Regression and adversarial checks for P250, with no floating-point sign decisions."""
import json
import hashlib
from fractions import Fraction as F
from pathlib import Path
import unittest

from .adaptive_signed_tail import (BETA_CUTOFF, CERT_CUTOFF, NEGATIVE_WITNESS,
                                   SignedSymbol, adaptive_search, load_p245b, run)

EXPECTED_CERT_SHA256='2bde93ccce25b8d2a1c6e59b49bcce43475fb2c7bc83366edaa1a19b7d4efe38'

class P250Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.arith=load_p245b()
        cls.symbol=SignedSymbol(cls.arith)
        cls.cert=run()

    def test_full_certification_replays_deterministically(self):
        self.assertEqual(run(), self.cert)
        wire=json.dumps(self.cert,indent=2,ensure_ascii=False)+"\n"
        self.assertEqual(hashlib.sha256(wire.encode()).hexdigest(), EXPECTED_CERT_SHA256)

    def test_exact_cover_no_missing_or_overlapping_cells(self):
        leaves=self.cert['leaves']
        self.assertEqual(F(leaves[0][0]),F(CERT_CUTOFF))
        self.assertEqual(F(leaves[-1][1]),F(BETA_CUTOFF))
        self.assertEqual(len(leaves),self.cert['count_cells'])
        for i,(a,b,z) in enumerate(leaves):
            self.assertLess(F(a),F(b))
            self.assertGreater(int(z),0)
            if i:
                self.assertEqual(F(a),F(leaves[i-1][1]))

    def test_bad_coverage_is_rejected_by_reconstruction(self):
        leaves=self.cert['leaves'][:]
        leaves.pop(2)
        self.assertNotEqual(F(leaves[2][0]),F(leaves[1][1]))

    def test_signed_symbol_witness_is_strictly_negative(self):
        signed=self.symbol.true_symbol_at(NEGATIVE_WITNESS)
        self.assertLess(signed.u,0)
        self.assertEqual(signed.desc(),self.cert['negative_actual_symbol_interval'])

    def test_positive_at_309_and_binet_cutoff(self):
        self.assertGreater(self.symbol.true_symbol_at(F(309)).l,0)
        self.assertEqual(self.cert['certified_positive_frequency_range'],[309,'infinity'])

    def test_coarse_cell_requires_refinement(self):
        self.assertLessEqual(self.symbol.gap_lower(F(309),F(504)).l,0)
        self.assertEqual(adaptive_search(self.symbol)['cells'],230)

    def test_invalid_analytic_range_fails_closed(self):
        with self.assertRaises(ValueError):
            self.symbol.gap_lower(F(308),F(505))

    def test_exact_negative_at_30845_rules_out_all_integer_cuts_below_309(self):
        self.assertTrue(F(308)<NEGATIVE_WITNESS<F(309))
        self.assertLess(self.symbol.true_symbol_at(NEGATIVE_WITNESS).u,0)
        self.assertEqual(self.cert['smallest_integer_tail_cutoff_for_this_fixed_beta'],309)

if __name__=='__main__': unittest.main()
