#!/usr/bin/env python3
"""Independent structural checks, run normally and with Python -O."""
import unittest
import json
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
from verify_spectral_rational import (read_zero_intervals,SCALING,cos_iv,
    pi_iv,exp_up_rational,positive_height_tail_bound,certify)

FILE=Path(__file__).with_name('preview_zeros.json')
DATA=json.loads(FILE.read_text())

class TestP72(unittest.TestCase):
    def test_zero_table_count_and_order(self):
        intervals=read_zero_intervals(DATA,20)
        self.assertEqual(len(intervals),20)
        self.assertTrue(all(g[0]<g[1] for g in intervals))

    def test_cosine_outward_inclusion(self):
        mp.mp.dps=65
        import mpmath
        for x in ('0','-3.9','1.1','3.1'):
            from fractions import Fraction as F
            c=cos_iv((int(F(x)*SCALING),int(F(x)*SCALING)))
            actual=mp.cos(mp.mpf(x))
            self.assertLessEqual(mp.mpf(c[0])/SCALING,actual)
            self.assertGreaterEqual(mp.mpf(c[1])/SCALING,actual)

    def test_exponential_upper(self):
        mp.mp.dps=65
        for v in ('20','22','23.74'):
            q=exp_up_rational(Q(v))
            self.assertGreater(mp.mpf(q.numerator)/q.denominator,mp.exp(mp.mpf(v)))

    def test_zero_tail_analytic(self):
        self.assertLess(positive_height_tail_bound(),Q(6,3*10**12))

    def test_example_spectral_point_against_mpmath(self):
        mp.mp.dps=55
        for n in (3,10,20):
            rec=certify(DATA,n,18005,5)
            self.assertTrue(rec['positive_given_input_intervals_are_valid'])
            self.assertFalse(rec['certified_positive'])
            low=Q(int(rec['grid_min_scaled']),SCALING)
            vals=[mp.im(mp.zetazero(i)) for i in range(1,n+1)]
            direct=min(sum(2*(1-mp.cos(g*t))/g**2 for g in vals) for t in (mp.mpf(18),mp.mpf('18.005')))
            self.assertLess(mp.mpf(low.numerator)/low.denominator,direct)

if __name__=='__main__':unittest.main()
