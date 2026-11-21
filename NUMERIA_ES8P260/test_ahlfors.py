"""Exact-rational regression and adversarial tests for NUMERIA ES8-P260."""
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest
from . import ahlfors_digamma as a

class AhlforsDigammaTests(unittest.TestCase):
    def test_prime_powers_exact(self):
        self.assertEqual(a.PRIME_POWERS,((2,2),(3,3),(4,2),(5,5)))
    def test_pi_bracket_strict(self):
        lo,hi=a.pi_interval()
        self.assertTrue(F(314159,100000)<lo<hi<F(355,113))
    def test_log_interval_simple_brackets(self):
        for x in (F(2),F(3),F(5),F(355,113)):
            lo,hi=a.log_interval(x)
            self.assertLess(lo,hi)
            self.assertLess(hi-lo,F(1,10**20))
    def test_log_invalid(self):
        for x in (1,F(1,2),-1):
            with self.assertRaises(ValueError):a.log_interval(x)
    def test_sqrt_grid(self):
        for q in (2,3,4,5):
            s=a.sqrt_lower(q)
            self.assertLessEqual(s*s,q)
            self.assertGreater((s+F(1,a.SCALE))**2,q)
    def test_sqrt_invalid(self):
        for q in (0,-1,2.0):
            with self.assertRaises(ValueError):a.sqrt_lower(q)
    def test_gamma_upper(self):
        g=a.euler_gamma_upper()
        self.assertTrue(F(5772,10000)<g<F(5774,10000))
    def test_harmonic_input_invalid(self):
        for n in (0,3,-1):
            with self.assertRaises(ValueError):a.euler_gamma_upper(n)
    def test_archimedean_anchor_strict(self):
        anchor=a.archimedean_anchor_lower()
        self.assertGreater(anchor,F(1,400))
        self.assertLess(anchor,F(1,200))
    def test_monotone_increment_lower(self):
        ds=[a.digamma_increment_lower(t) for t,_ in a.FLOORS]
        self.assertEqual(len(ds),4)
        self.assertTrue(F(0)<ds[0]<ds[1]<ds[2]<ds[3])
    def test_floor_850(self):
        self.assertGreater(a.archimedean_anchor_lower()+a.digamma_increment_lower(850),F(1,2))
    def test_new_floor_1100(self):
        self.assertGreater(a.archimedean_anchor_lower()+a.digamma_increment_lower(1100),F(3,4))
    def test_new_floor_1500(self):
        self.assertGreater(a.archimedean_anchor_lower()+a.digamma_increment_lower(1500),F(1))
    def test_new_floor_2000(self):
        self.assertGreater(a.archimedean_anchor_lower()+a.digamma_increment_lower(2000),F(5,4))
    def test_false_higher_floor_rejected(self):
        lower=a.archimedean_anchor_lower()+a.digamma_increment_lower(850)
        with self.assertRaises(ArithmeticError):
            a.require(lower>F(3,4),'invalid stricter floor at 850')
    def test_delta_wrong_input_rejected(self):
        for t in (503,504,505,850.0):
            with self.assertRaises(ValueError):a.digamma_increment_lower(t)
    def test_json_honest_flags(self):
        j=a.audit()
        self.assertFalse(j['sector_odd_377_positive_certified'])
        self.assertFalse(j['Weil_crossing_established'])
        self.assertFalse(j['Riemann_hypothesis_proved'])
        self.assertEqual(len(j['new_piecewise_floors']),4)
        self.assertEqual(j['new_piecewise_floors'][1]['certified_constant_floor_exact'],'3/4')
    def test_underflow_cannot_flip_sign(self):
        for x in (F(-1,10**40),F(0),F(1,10**40)):
            self.assertLessEqual(a.conservative_grid_lower(x),x)
    def test_cli_replay(self):
        with TemporaryDirectory() as folder:
            out=Path(folder)/'r.json'
            out.write_text(json.dumps(a.audit(),indent=2,ensure_ascii=False)+'\n')
            obj=json.loads(out.read_text())
            self.assertEqual(obj['schema'],'NUMERIA.ES8P260.AHLFORS_ARCHIMEDEAN_RATIONAL.1')

if __name__=='__main__':unittest.main()
