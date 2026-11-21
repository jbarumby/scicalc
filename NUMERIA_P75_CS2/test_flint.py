"""Arb-only regression gates: run in GitHub Actions with python-flint==0.9.0."""
from fractions import Fraction as F
import unittest
from flint import arb, ctx
from NUMERIA_P75_CS2.cs_v2 import TABLE_V2, row_indices, density_integral_upper, ChourasiyaV2Tail
from NUMERIA_P75_CS2.local import local_curvature_upper, lower_secant_local
from NUMERIA_P74_ZERO_SIDE.zeros import verified_ordinates
from NUMERIA_P74_ZERO_SIDE.certify import Spectrum, ball
from NUMERIA_P75_CS2.certify import run

class TestDirectedArb(unittest.TestCase):
    def test_source_table(self):
        self.assertEqual(len(TABLE_V2),16)
        self.assertEqual(TABLE_V2[0],('5.069','9.461','167.8'))
        self.assertEqual(TABLE_V2[-1],('46.06','4.897','86.84'))
        self.assertEqual(row_indices(F(17,32)),(0,1))
        self.assertEqual(row_indices(F(1)),(15,))

    def test_density_monotonic_caps(self):
        with ctx.workprec(192):
            tail=ChourasiyaV2Tail(32)
            self.assertEqual(len(tail.mass_caps),32)
            self.assertTrue(all(a>=b for a,b in zip(tail.mass_caps,tail.mass_caps[1:])))
            self.assertTrue(tail.bound(arb(79))>0)
            self.assertTrue(density_integral_upper(F(3,4))>0)

    def test_refinement(self):
        with ctx.workprec(192):
            coarse=ChourasiyaV2Tail(16)
            fine=ChourasiyaV2Tail(32)
            for t in (arb(60),arb(79)):
                # Tiny differences due to outward rounding are allowed.
                self.assertTrue(fine.bound(t)<coarse.bound(t)+arb('1e-20'))

    def test_local_bound_reduces_to_global(self):
        with ctx.workprec(192):
            zeros,meta=verified_ordinates(3)
            width=arb(1)/64
            maximum=local_curvature_upper(zeros,arb(18)+width/2,width)
            self.assertTrue(maximum<=arb(6))
            self.assertTrue(maximum>0)
            lower=lower_secant_local(arb(1),arb(1),width,maximum)
            global_lower=arb(1)-arb(3)*width*width/4
            self.assertTrue(lower>=global_lower-arb('1e-40'))

    def test_full_short_interval(self):
        r=run(end='18.03125',count=3,prec=192,segments=16,mesh='1/64')
        self.assertEqual(r['schema'],'NUMERIA.P75.CS_V2.LOCAL_CURVATURE.1')
        self.assertTrue(r['overlaps_p71'])
        self.assertEqual(r['tails']['segments'],16)
        self.assertIn(r['status'],('conditional_finite_positive','inconclusive'))

if __name__=='__main__':unittest.main()
