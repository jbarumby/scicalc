#!/usr/bin/env python3
"""Regression tests; no assertions used for gates in production scripts."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import unittest

from certify_adaptive import certify, spectral_lower
from experiment import discover, ordinates

ROOT=Path(__file__).resolve().parents[1]
ZEROS=ROOT/'NUMERIA_P72_ZERO_SIDE'/'arb_zeros_20.json'
EXPECTED_HASH='c5993b75a6db7803251655955859cd935f4381b2c0bc20fd77b9f0258bdeddab'

class P73Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=ZEROS.read_bytes()
        cls.data=json.loads(raw)
        cls.sha=sha256(raw).hexdigest()

    def test_input_hash_pin(self):
        self.assertEqual(self.sha,EXPECTED_HASH)
        self.assertEqual(len(self.data['zeros']),20)

    def test_cosine_at_zero_positive(self):
        import verify_spectral_rational as p72
        ordinates=p72.read_zero_intervals(self.data,3)
        lo=spectral_lower(F(0),ordinates,p72.pi_iv())
        self.assertLessEqual(lo,F(0))
        self.assertGreater(lo,-F(1,10**40))

    def test_rational_certification_is_reproducible(self):
        a=certify(self.data,3,F(18),F('18.2'),bucket=F('0.05'))
        b=certify(self.data,3,F(18),F('18.2'),bucket=F('0.05'))
        self.assertEqual(a,b)
        self.assertTrue(a['strictly_positive_conditional_on_zero_data_and_external_theorems'])
        self.assertGreater(F(*map(int,a['min_lower_rational'])),0)

    def test_20_zeros_to_46_6_is_positive_and_replay_matches(self):
        expected=json.loads((Path(__file__).with_name('certificate_20_t46.6.json')).read_text())
        result=certify(self.data,20,F(18),F('46.6'),F('0.05'))
        result['zeros_file_sha256']=self.sha
        self.assertEqual(result,expected)
        self.assertGreater(F(*map(int,result['min_lower_rational'])),0)
        self.assertEqual(result['cell_count'],565)

    def test_negative_tail_proxy_is_not_zeta_counterexample(self):
        gammas=ordinates(self.data,20)
        a=discover('uniform',gammas,budget=128)
        self.assertGreater(a['negative_samples'],0)
        self.assertLess(a['minimum_heuristic_score'],0)

    def test_random_reproducibility(self):
        g=ordinates(self.data,20)
        self.assertEqual(discover('random',g,budget=128,seed=73),
                         discover('random',g,budget=128,seed=73))

    def test_reject_tampered_input_sha(self):
        runner=Path(__file__).with_name('certify_adaptive.py')
        run=subprocess.run([sys.executable,str(runner),'--zeros',str(ZEROS),
                  '--expected-zeros-sha256','0'*64,'--to','18.1'],
                  capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('SHA mismatch',run.stderr+run.stdout)

if __name__=='__main__':unittest.main()
