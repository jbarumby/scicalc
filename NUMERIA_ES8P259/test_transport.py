"""Exact rational external postcheck of GitHub Arb-generated P259 certificate."""
import json
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from .transport import (audit, log7_bounds, verify_derivative_budget,
                        replay_result, L_TEST, L_BELOW,
                        EXPECTED_GATE_B_SHA256, NEW_DILATION_LIPSCHITZ)

SOURCE=Path('replay/P259_GATE_B_PREFIX.json')

class GateBTransport(unittest.TestCase):
    def test_two_sided_threshold(self):
        lo,hi=log7_bounds()
        self.assertTrue(L_BELOW<lo/2<hi/2<L_TEST)

    def test_uniform_derivative_budget(self):
        x=verify_derivative_budget()
        self.assertEqual(x['P258_existing_band_and_poles_upper_exact'],'21246')
        self.assertEqual(x['P259_new_850_highpass_derivative_upper_exact'],'851/2')
        self.assertEqual(NEW_DILATION_LIPSCHITZ,21676)

    def test_pinned_source(self):
        self.assertEqual(len(EXPECTED_GATE_B_SHA256),64)

    def test_certified_original_and_transport(self):
        if not SOURCE.is_file():
            self.skipTest('P259 Gate B artifact has not been retrieved')
        a=audit(SOURCE)
        self.assertTrue(a['same_fixed_piecewise_four_prime_auxiliary_not_positive_for_any_L_ge_log7_half'])
        self.assertLess(F(a['transported_upper_exact']),0)
        self.assertLess(F(a['loss_exact']),F(1,10**27))
        self.assertLess(F(a['exact_original_rayleigh_upper']),F(-4,10**27))

    def test_tampered_source_rejected(self):
        if not SOURCE.is_file():
            self.skipTest('P259 Gate B artifact has not been retrieved')
        with TemporaryDirectory() as d:
            bad=Path(d)/'tampered.json'
            bad.write_bytes(SOURCE.read_bytes()+b'bad')
            with self.assertRaisesRegex(ArithmeticError,'SHA256'):
                audit(bad)

    def test_metadata_mutation_rejected(self):
        if not SOURCE.is_file():
            self.skipTest('P259 Gate B artifact has not been retrieved')
        j=json.loads(SOURCE.read_text())
        j['checked_prefixes'][-1]['negative_witness']['exact_rational_replay']['operator_norm_errors_exact']['P253_weighted_GL60']='0'
        with self.assertRaises(ArithmeticError):
            replay_result(j)

    def test_positive_forgery_rejected(self):
        if not SOURCE.is_file():
            self.skipTest('P259 Gate B artifact has not been retrieved')
        j=json.loads(SOURCE.read_text())
        j['status']='CERTIFIED_ODD_HEAD_POSITIVE'
        with self.assertRaises(ArithmeticError):
            replay_result(j)

if __name__=='__main__':unittest.main()
