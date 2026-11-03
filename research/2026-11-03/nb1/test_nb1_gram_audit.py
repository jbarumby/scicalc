#!/usr/bin/env python3
from fractions import Fraction as F
import json
import pathlib
import unittest
import mpmath as mp
from nb1_gram_audit import cuts, certify, log_rational_enclosure, E

class GramCertificateTests(unittest.TestCase):
    def test_zero_coefficients_recover_unit_norm_upper(self):
        c=certify(1,cuts(1),[F(0)])
        self.assertGreaterEqual(F(c['squared_error_upper_bound']), F(1))
        self.assertLess(F(c['squared_error_upper_bound']), F(1)+F(100000,10**26))
        self.assertLessEqual(F(c['squared_error_lower_bound']), F(1)-E)

    def test_piecewise_representation(self):
        for N in (6,8,10):
            for k in (1,2,7,13,105,1000,9999):
                x=(F(1,k+1)+F(1,k))/2
                for n in range(1,N+1):
                    z=F(1,n)/x
                    rho=z-(z.numerator//z.denominator)
                    self.assertEqual(rho,F(1,n)/x-(k//n))

    def test_rational_log_brackets(self):
        mp.mp.dps=90
        for k in (1,2,5,50,9999):
            r=F(k+1,k)
            lo,hi=log_rational_enclosure(r)
            val=mp.log(mp.mpf(k+1)/k)
            self.assertLessEqual(mp.mpf(lo.numerator)/lo.denominator,val)
            self.assertLessEqual(val,mp.mpf(hi.numerator)/hi.denominator)

    def test_zero_padding_invariance(self):
        record=json.loads(pathlib.Path(__file__).with_name('RESULTADOS_NB1.json').read_text())['n_results']
        c=[F(q) for q in record[1]['coefficients']]
        a=certify(8,cuts(8),c)
        b=certify(10,cuts(10),c+[F(0),F(0)])
        self.assertEqual(a['squared_error_upper_bound'],b['squared_error_upper_bound'])
        self.assertEqual(a['squared_error_lower_bound'],b['squared_error_lower_bound'])

    def test_recorded_ranges(self):
        rows=json.loads(pathlib.Path(__file__).with_name('RESULTADOS_NB1.json').read_text())['n_results']
        self.assertEqual([r['n'] for r in rows],[6,8,10])
        for r in rows:
            self.assertLess(F(0),F(r['squared_error_lower_bound']))
            self.assertLess(F(r['squared_error_lower_bound']),F(r['squared_error_upper_bound']))
            self.assertLess(F(r['squared_error_upper_bound']),F(1))

if __name__ == '__main__': unittest.main(verbosity=2)
