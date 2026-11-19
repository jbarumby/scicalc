#!/usr/bin/env python3
"""Exact and adversarial tests for the P245 Gate B interval quadrature pilot."""
import copy
import hashlib
import json
import math
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path

from p245b_certified_gl60 import (I,PI,PREC,S,ROOT_PREC,ROOT_SCALE,
    bernoulli_even,certified_nodes,legendre_at,log_interval,iv,re_digamma,
    sin_cos,spherical_j_small,bessel_small_arg_series,verify_against_pilot,
    source_ball_to_fraction,PILOT_HASHES,PAIRS)

HERE=Path(__file__).resolve().parent
ROOTS=HERE/'P245B_ROOT_BRACKETS.json'
PILOTS=HERE/'pilot_sources'
REPORT=HERE/'P245B_CERT.json'

class TestP245BGates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(REPORT.read_text())
        cls.rootdata=json.loads(ROOTS.read_text())

    def test_01_60_root_intervals_by_exact_sign(self):
        roots=certified_nodes(ROOTS)
        self.assertEqual(len(roots),60)
        for i,r in enumerate(self.rootdata['roots']):
            a=Q(int(r['lower_numerator']),ROOT_SCALE)
            b=Q(int(r['upper_numerator']),ROOT_SCALE)
            self.assertEqual(b-a,Q(1,ROOT_SCALE))
            self.assertLess(legendre_at(a,60)*legendre_at(b,60),0)
            self.assertGreater(roots[i][1].l,0)
        self.assertTrue(sum((w for x,w,m in roots),iv(0)).contain_fraction(Q(2)))

    def test_02_const_pi_bound(self):
        self.assertLess(3*S,PI.l)
        self.assertLess(PI.u*7,22*S)
        self.assertTrue(log_interval(2).contain_interval(log_interval(iv(2))))

    def test_03_bernoulli_index_regression(self):
        b=bernoulli_even()
        self.assertEqual(b[0],Q(1,6))
        self.assertEqual(b[1],Q(-1,30))
        self.assertEqual(b[3],Q(-1,30))
        self.assertEqual(b[2],Q(1,42))
        self.assertNotEqual(Q(-b[1],4),Q(-b[3],8))

    def test_04_trig_known_values(self):
        sn,cs=sin_cos(PI)
        self.assertTrue(sn.contain_fraction(Q(0)))
        self.assertTrue(cs.contain_fraction(Q(-1)))
        sn,cs=sin_cos(PI/2)
        self.assertTrue(sn.contain_fraction(Q(1)))
        self.assertTrue(cs.contain_fraction(Q(0)))

    def test_05_spherical_bessel_low_degree_agrees_with_series(self):
        x=iv(Q(1,2))
        v=spherical_j_small(x)
        for n in (2,3):
            s=bessel_small_arg_series(x,n)
            self.assertLessEqual(max(v[n].l,s.l),min(v[n].u,s.u))

    def test_06_digamma_interval_finite(self):
        for y in (iv(0),iv(1),iv(252)):
            a=re_digamma(y)
            self.assertLess(a.l,a.u)
            self.assertLess(abs(a.l),10*S)
            self.assertLess(abs(a.u),10*S)

    def test_07_source_sha256_pins(self):
        self.assertEqual(len(PILOT_HASHES),6)
        for f,h in PILOT_HASHES.items():
            self.assertEqual(hashlib.sha256((PILOTS/f).read_bytes()).hexdigest(),h)

    def test_08_all_48_pilot_enclosures_and_groups(self):
        self.assertEqual(self.data['certified_selected_entries'],48)
        self.assertEqual(len(self.data['pilots']),6)
        self.assertEqual(self.data['parameters']['precision_bits'],256)
        for row in self.data['pilots']:
            self.assertEqual(len(row['comparisons']),8)
            self.assertEqual(len(row['high_degree_magnitude_bounds']),3)
            # Re-evaluate source containment using exact Fraction, not booleans.
            source=json.loads((PILOTS/row['source_file']).read_text())
            samples={(i,j):b for i,j,b in source['gauss_interval_matrix_sample']}
            for sample in row['comparisons']:
                lo,hi=source_ball_to_fraction(samples[tuple(sample['pair'])])
                b=sample['certified_interval']
                v=I(int(b['lower']),int(b['upper']))
                self.assertTrue(v.contain_fraction(lo))
                self.assertTrue(v.contain_fraction(hi))

    def test_09_large_degree_rational_bounds(self):
        for row in self.data['pilots']:
            m={tuple(e['pair']):e['absolute_value_strictly_below_10_to_negative_power']
               for e in row['high_degree_magnitude_bounds']}
            if row['panel']==1007:
                self.assertGreaterEqual(m[0,376],102)
                self.assertGreaterEqual(m[376,376],206)
            if row['panel']==0:
                self.assertGreaterEqual(m[0,376],2360)
            for e in row['high_degree_magnitude_bounds']:
                self.assertTrue(e['verified_exact_rational'])

    def test_10_tampered_root_rejected(self):
        x=copy.deepcopy(self.rootdata)
        x['roots'][0]['lower_numerator']=str(int(x['roots'][0]['lower_numerator'])+3)
        x['roots'][0]['upper_numerator']=str(int(x['roots'][0]['upper_numerator'])+3)
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/'roots.json';file.write_text(json.dumps(x))
            with self.assertRaises(ArithmeticError):certified_nodes(file)

    def test_11_tampered_pilot_rejected_by_hash(self):
        row=self.data['pilots'][0]
        with tempfile.TemporaryDirectory() as td:
            td=Path(td)
            source=PILOTS/row['source_file']
            d=json.loads(source.read_text())
            d['gauss_interval_matrix_sample'][0][2]=[['0','0'],['0','0']]
            (td/source.name).write_text(json.dumps(d))
            values={tuple(z['pair']):I(int(z['certified_interval']['lower']),
                                      int(z['certified_interval']['upper']))
                    for z in row['comparisons']}
            with self.assertRaises(ArithmeticError):
                verify_against_pilot(row['panel'],int(row['sector']=='odd'),values,td)

    def test_12_status_no_global_proof(self):
        self.assertFalse(self.data['all_GL60_panels_certified'])
        self.assertFalse(self.data['Weil_Zhu_identity_certified'])
        self.assertFalse(self.data['RH_proved'])

if __name__=='__main__':unittest.main()
