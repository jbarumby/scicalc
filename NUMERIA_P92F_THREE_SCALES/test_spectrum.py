"""Exact verification of refined scale-spectrum Schur envelopes."""
import unittest
from fractions import Fraction as F
from .spectrum_gate import certify,spectrum_rows,TARGET_LOWER
from .certificate import points_deltas_from_gaps

class SpectrumTests(unittest.TestCase):
    def test_exact_five_square_certificate(self):
        c=certify()
        self.assertEqual(c['five_scale_rational_norm_upper'],'6829/2160')
        self.assertTrue(c['five_scales_c_times_1_4_9_16_25_unlimited_points_proved'])
        self.assertLess(F(6829,2160),TARGET_LOWER)
        self.assertFalse(c['all_arbitrary_scale_spectra_global_Delta_proved'])
        self.assertFalse(c['Riemann_Hypothesis_proved'])

    def test_at_least_one_configuration_realizes_all_five_scales(self):
        gaps=(1,4,4,9,9,16,16,25,25)
        _,d=points_deltas_from_gaps(gaps)
        self.assertEqual(set(d),set(map(F,(1,4,9,16,25))))

    def test_rational_scaling_invariance(self):
        a=max(F(row['rational_upper']) for row in spectrum_rows((1,2,3,4,5)))
        b=max(F(row['rational_upper']) for row in spectrum_rows(tuple(F(t)*F(3,2) for t in (1,2,3,4,5))))
        self.assertEqual(a,b)

    def test_the_test_rejects_far_apart_four_scales(self):
        c=certify()
        bad=next(x for x in c['cases'] if x['name']=='four_extreme')
        self.assertFalse(bad['P88_target_proved_by_rational_lower'])
        self.assertGreater(F(bad['max_rational_schur_ceiling']),TARGET_LOWER)

    def test_close_six_scales_are_certifiable(self):
        c=certify()
        good=next(x for x in c['cases'] if x['name']=='six_close')
        self.assertTrue(good['P88_target_proved_by_rational_lower'])

    def test_invalid_spectra(self):
        for values in ((),(0,),(-1,2)):
            with self.assertRaises(ValueError):spectrum_rows(values)

if __name__=='__main__':
    unittest.main()
