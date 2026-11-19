#!/usr/bin/env python3
"""Adversarial exact-arithmetic checks for P245 Gate C."""
import copy, gzip, hashlib, json, tempfile, unittest
from pathlib import Path
from fractions import Fraction as Q
from p245b_certified_gl60 import I, S, certified_nodes, source_ball_to_fraction
from p245c_chunk25 import HASHES, PAIRS, N, triangular_pos, bessel_series, get_source
from p245c_gl60_error_budget import run as budget_run, amplitudes, ceil_sqrt_q

HERE=Path(__file__).resolve().parent

class TestP245C(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=json.loads((HERE/'P245C_CERT.json').read_text())
        cls.error=json.loads((HERE/'P245C_GL60_ERROR_BUDGET.json').read_text())

    def test_01_source_sha256_and_25_panel_metadata(self):
        for sec in HASHES:
            d,p=get_source(HERE,sec)
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),HASHES[sec])
            self.assertEqual(d['GL60_node_evaluations'],1500)
            self.assertEqual(d['upper_triangle_count'],N*(N+1)//2)

    def test_02_60_roots_reproved_without_approximate_root_finder(self):
        roots=certified_nodes(HERE/'P245B_ROOT_BRACKETS.json')
        self.assertEqual(len(roots),60)
        self.assertGreater(min(w.l for _,w,_ in roots),0)

    def test_03_all_26_chunk_sum_containments_exact(self):
        self.assertEqual(self.cert['verified_chunk_sum_enclosures'],2*len(PAIRS))
        self.assertEqual(self.cert['all_panel_pair_evaluations'],25*2*len(PAIRS))
        for rec in self.cert['chunk_sums']:
            d,_=get_source(HERE,rec['sector'])
            self.assertEqual(len(rec['interval_pairs']),len(PAIRS))
            for item in rec['interval_pairs']:
                i,j=item['pair']
                ii,jj,ball=d['C_partial_upper_triangle'][triangular_pos(i,j)]
                self.assertEqual((ii,jj),(i,j))
                lo,hi=source_ball_to_fraction(ball)
                saved=item['interval']
                enclosing=I(int(saved['lower']),int(saved['upper']))
                self.assertTrue(enclosing.contain_fraction(lo))
                self.assertTrue(enclosing.contain_fraction(hi))

    def test_04_error_budget_recomputed_exact(self):
        renewed=budget_run()
        self.assertEqual(renewed,self.error)
        for sec in renewed['sectors'].values():
            self.assertTrue(sec['max_row_error_less_than_4e-47'])

    def test_05_754_degree_amplitudes_and_bounds(self):
        for p in (0,1):
            a=amplitudes(p)
            self.assertEqual(len(a),377)
            self.assertTrue(all(0<x for x in a))
            self.assertLessEqual(max(a),100)

    def test_06_roots_tampered_are_rejected(self):
        x=json.loads((HERE/'P245B_ROOT_BRACKETS.json').read_text())
        x['roots'][0]['lower_numerator']=str(int(x['roots'][0]['lower_numerator'])+5)
        x['roots'][0]['upper_numerator']=str(int(x['roots'][0]['upper_numerator'])+5)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'P245B_ROOT_BRACKETS.json';p.write_text(json.dumps(x))
            with self.assertRaises(ArithmeticError): certified_nodes(p)

    def test_07_source_mutation_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td)
            d=(HERE/'P236_CHUNK25_even.json.gz').read_bytes()
            (td/'P236_CHUNK25_even.json.gz').write_bytes(d+b'\0')
            with self.assertRaises(ArithmeticError):get_source(td,'even')

    def test_08_forged_report_interval_does_not_contain_source(self):
        item=self.cert['chunk_sums'][0]['interval_pairs'][0]
        d,_=get_source(HERE,'even')
        i,j=item['pair']; _,_,ball=d['C_partial_upper_triangle'][triangular_pos(i,j)]
        low,high=source_ball_to_fraction(ball)
        fake=I(0,0)
        self.assertFalse(fake.contain_fraction(low) and fake.contain_fraction(high))

    def test_09_bessel_series_consistent_with_elementary_low_degrees(self):
        from p245b_certified_gl60 import spherical_j_small,iv
        for x in [Q(1,2),Q(5)]:
            ref=spherical_j_small(iv(x))
            for n in range(4):
                val=bessel_series(iv(x),n)
                self.assertLessEqual(max(val.l,ref[n].l),min(val.u,ref[n].u))

    def test_10_square_root_ceil_exact(self):
        self.assertEqual(ceil_sqrt_q(Q(9,4)),2)
        self.assertEqual(ceil_sqrt_q(Q(4)),2)
        self.assertEqual(ceil_sqrt_q(Q(101,100)),2)

    def test_11_scope_restrictions_preserved(self):
        self.assertFalse(self.cert['actual_integrals_fully_certified'])
        self.assertFalse(self.cert['all_142506_entries_certified_independently'])
        self.assertFalse(self.cert['RH_proved'])
        self.assertFalse(self.error['complex_analytic_hypotheses_independently_proved'])

    def test_12_stratified_156_entries_and_external_overlap(self):
        cert=json.loads((HERE/'P245C_STRATIFIED_CERT.json').read_text())
        self.assertEqual(cert['panel_entries_reconstructed'],156)
        self.assertEqual(cert['externally_compared_panel_entries'],52)
        self.assertEqual(len(cert['panels']),12)
        self.assertEqual(sorted({r['panel'] for r in cert['panels']}),[100,250,503,750,875,1007])
        for block in cert['panels']:
            self.assertEqual(len(block['entries']),len(PAIRS))
            for item in block['entries']:
                if item['external_Arb_interval_compared']:
                    self.assertTrue(item['external_Arb_interval_contained'])
            if block['panel'] not in (503,1007):
                self.assertIsNone(block['P236_source'])

    def test_13_stratified_external_source_interval_containment_exact(self):
        cert=json.loads((HERE/'P245C_STRATIFIED_CERT.json').read_text())
        from p245b_certified_gl60 import PILOT_HASHES
        comparisons=0
        for block in cert['panels']:
            if block['P236_source'] is None:continue
            name=block['P236_source']
            f=HERE/'pilot_sources'/name
            self.assertEqual(hashlib.sha256(f.read_bytes()).hexdigest(),PILOT_HASHES[name])
            src=json.loads(f.read_text())
            source={(i,j):ball for i,j,ball in src['gauss_interval_matrix_sample']}
            for entry in block['entries']:
                i,j=entry['pair']
                if (i,j) not in source:continue
                lo,hi=source_ball_to_fraction(source[i,j])
                b=entry['interval'];ival=I(int(b['lower']),int(b['upper']))
                self.assertTrue(ival.contain_fraction(lo) and ival.contain_fraction(hi))
                comparisons+=1
        self.assertEqual(comparisons,52)

    def test_14_pin_the_b_gate_source_in_new_archive(self):
        p=HERE/'p245b_certified_gl60.py'
        self.assertGreater(p.stat().st_size,10000)
        self.assertTrue((HERE/'P245B_ROOT_BRACKETS.json').exists())

if __name__=='__main__':unittest.main()
