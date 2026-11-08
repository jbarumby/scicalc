#!/usr/bin/env python3
"""ES8-P102: inexpensive FLINT 0.9.0 interval callback contract tests.
Passing these tests is NOT a complete audit of nested acb.integral.
"""
import sys, unittest
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'es8_p99'))
from flint import arb, acb, fmpq, ctx
from W12_arb_direto_P98 import W12P92
ctx.prec=256
def frac(v):
    q=v.fmpq()
    return F(int(q.p),int(q.q))
class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=W12P92.make('1e-8',1000,12,eps_denominator=2048)
    def test_rational_interval_enclosure(self):
        cases=[(F(1,3),F(2,3)),(F(-1,7),F(1,11)),
               (F(0),F(0)),(F(-11,20),F(-9,20)),
               (F(1,2**201),F(1,2**200))]
        for lo,hi in cases:
            z=self.c.ball_from_rational_interval(lo,hi)
            self.assertLessEqual(frac(z.lower()),lo)
            self.assertGreaterEqual(frac(z.upper()),hi)
    def test_bump_crossing_support_not_analytic(self):
        h=arb(fmpq(2,5))
        self.assertFalse(self.c.bump(acb(arb(h,arb('0.01'))),True).is_finite())
    def test_real_bump_range_encloses_zero_and_max(self):
        b=self.c.bump(acb(arb(0,arb(fmpq(1,2)))),False)
        self.assertTrue(b.is_finite())
        self.assertLessEqual(frac(b.real.lower()),0)
        self.assertGreaterEqual(frac(b.real.upper()),F(1,4))
    def test_bump_inside_analytic(self):
        self.assertTrue(self.c.bump(acb(arb(0,arb(fmpq(1,10)))),True).is_finite())
    def test_zero_beyond_support(self):
        self.assertTrue(self.c.bump(acb(arb(2)),True).is_zero())
        self.assertTrue(self.c.correlation(acb(arb(2))).is_zero())
    def test_regularized_kernel(self):
        for t in (arb(0),arb(fmpq(1,100)),arb(fmpq(1,2)),arb(1)):
            y=(-t/2).exp()
            stable=2*y/((1+y)*(1+y*y))
            if t.is_zero():self.assertTrue(stable.contains(arb(fmpq(1,2))))
            else:self.assertTrue(stable.overlaps(((t/2).exp()-1)/t.sinh()))
    def test_exact_origin_bound(self):
        e=F(1,2048);h=F(2,5)
        self.assertEqual(e*e/h,F(5,8388608))
    def test_strict_derivative_bound(self):
        self.assertTrue(bool(16/(arb.const_e()**3)<2))
    def test_empty_intersection_rejected(self):
        with self.assertRaises(ArithmeticError):
            self.c.ball_from_rational_interval(F(2),F(1))
if __name__=='__main__':unittest.main(verbosity=2)
