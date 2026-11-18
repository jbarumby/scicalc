#!/usr/bin/env python3
"""Lightweight P244 tests; the full MPFR replay is a separate explicit action."""
import json
import shutil
import subprocess
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path
from p244_analytic_gates import verify, exp_lower, exp_upper
from p244_replay_fail_closed import validate_replays, digest

ROOT=Path(__file__).resolve().parent
DATA=Path('/mnt/data')

class P244Tests(unittest.TestCase):
    def test_analytic_gates(self):
        report=verify()
        self.assertEqual(report['symbol_M20_rational_bound'],'1384/75')
        self.assertEqual(report['high_frequency_symbol_excess_strict_lower_bound'],'335/169344')
        self.assertFalse(report['RH_proved'])

    def test_exp_bounds(self):
        self.assertGreater(exp_lower(F(7,10),3),2)
        self.assertLess(exp_upper(F(161,100),8),6)
        self.assertGreater(exp_lower(F(161,100),8),5)

    def test_original_traces_bound(self):
        with tempfile.TemporaryDirectory() as d:
            for sec in ('even','odd'):
                n=f'P242_PIVOTS_{sec}_1024.txt'
                shutil.copy2(DATA/n,Path(d)/n)
            r=validate_replays(Path(d),DATA,DATA/'P242_AUDIT_CERT.json')
            self.assertEqual(len(r),2)
            self.assertTrue(r['even']['identical_reexecution'])

    def test_tampered_inverse_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            for sec in ('even','odd'):
                n=f'P242_PIVOTS_{sec}_1024.txt'
                raw=(DATA/n).read_text()
                if sec=='even':
                    raw='\n'.join('# inverse_frobenius_squared_upper 1 0' if ln.startswith('# inverse_frobenius_squared_upper ') else ln for ln in raw.splitlines())+'\n'
                (Path(d)/n).write_text(raw)
            with self.assertRaisesRegex(ArithmeticError,'regenerated trace SHA mismatch'):
                validate_replays(Path(d),DATA,DATA/'P242_AUDIT_CERT.json')

    def test_metadata_dimensions(self):
        r=json.loads((DATA/'P242_AUDIT_CERT.json').read_text())
        self.assertEqual(r['total_independent_positive_pivots'],754)
        self.assertFalse(r['full_Weil_Zhu_bridge_verified'])

    def test_numerical_spotcheck_is_diagnostic(self):
        r=json.loads((ROOT/'P244_NUMERICAL_SPOTCHECK.json').read_text())
        self.assertTrue(r['diagnostic_only_float64'])
        self.assertEqual(len(r['sectors']['odd']),6)
        self.assertEqual(len(r['sectors']['even']),6)
        for rows in r['sectors'].values():
            self.assertTrue(all(x['abs_diff']<1e-12 for x in rows))

    def test_C_MPFR_small_exact(self):
        # Independent exercise of the actual C algorithm on a rational 2x2.
        with tempfile.TemporaryDirectory() as d:
            exe=Path(d)/'ldl';p=Path(d)/'m';o=Path(d)/'p'
            c=subprocess.run(['gcc','-O2',str(ROOT/'p242_mpfr_interval_ldl.c'),'-o',str(exe),'-l:libmpfr.so.6','-lgmp'],capture_output=True,text=True)
            self.assertEqual(c.returncode,0,c.stderr)
            p.write_text('0 0 2 0\n0 1 1 0\n1 1 2 0\n')
            q=subprocess.run([str(exe),str(p),str(o),'256','2'],capture_output=True,text=True)
            self.assertEqual(q.returncode,0,q.stderr)
            rows=[x.split() for x in o.read_text().splitlines() if x and x[0].isdigit()]
            self.assertEqual(len(rows),2)
            for row,target in zip(rows,(F(2),F(3,2))):
                lo,hi=F(int(row[1]))*F(2)**int(row[2]),F(int(row[3]))*F(2)**int(row[4])
                self.assertLessEqual(lo,target);self.assertGreaterEqual(hi,target)

if __name__=='__main__':
    unittest.main(verbosity=2)
