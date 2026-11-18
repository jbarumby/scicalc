import unittest, subprocess, tempfile, sys
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from p242_verify_exact import dyad,ball_bounds,load_mpfr_trace

class TestNumeriaP242(unittest.TestCase):
    def test_dyadic_rational(self):
        self.assertEqual(dyad(-5,-3),Fraction(-5,8))
        self.assertEqual(dyad(2,3),16)
        self.assertEqual(ball_bounds([['3','-1'],['1','-2']]),(Fraction(5,4),Fraction(7,4)))

    @classmethod
    def setUpClass(cls):
        cls.exe=ROOT/'p242_mpfr_ldl'
        compile=['gcc','-O3','-Wall',str(ROOT/'p242_mpfr_interval_ldl.c'),'-o',str(cls.exe),'-l:libmpfr.so.6','-lgmp']
        r=subprocess.run(compile,capture_output=True,text=True)
        if r.returncode: raise RuntimeError(r.stderr)

    def test_exact_spd_2x2(self):
        # [[2,1],[1,2]]: LDL pivots 2,3/2, inverse F^2=10/9.
        with tempfile.TemporaryDirectory() as d:
            inp=Path(d)/'matrix.txt';out=Path(d)/'piv.txt'
            inp.write_text('0 0 2 0\n0 1 1 0\n1 1 2 0\n')
            r=subprocess.run([str(self.exe),str(inp),str(out),'256','2'],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            lines=[l for l in out.read_text().splitlines() if not l.startswith('#')]
            self.assertEqual(len(lines),2)
            for i,v in enumerate([Fraction(2),Fraction(3,2)]):
                p=lines[i].split()
                self.assertLessEqual(dyad(p[1],p[2]),v)
                self.assertGreaterEqual(dyad(p[3],p[4]),v)
            raw=out.read_text().split('# inverse_frobenius_squared_upper ')[1].splitlines()[0].split()
            nu=dyad(*raw)
            self.assertGreaterEqual(nu,Fraction(10,9))
            self.assertLess(nu,Fraction(112,100))

    def test_reject_nonpositive(self):
        with tempfile.TemporaryDirectory() as d:
            inp=Path(d)/'matrix.txt';out=Path(d)/'piv.txt'
            inp.write_text('0 0 -1 0\n0 1 0 0\n1 1 2 0\n')
            r=subprocess.run([str(self.exe),str(inp),str(out),'256','2'],capture_output=True,text=True)
            self.assertEqual(r.returncode,3)
            self.assertIn('INCONCLUSIVE pivot',r.stderr)

    def test_reject_bad_triangle(self):
        with tempfile.TemporaryDirectory() as d:
            inp=Path(d)/'matrix.txt';out=Path(d)/'piv.txt'
            inp.write_text('0 0 2 0\n1 0 1 0\n1 1 2 0\n')
            r=subprocess.run([str(self.exe),str(inp),str(out),'256','2'],capture_output=True,text=True)
            self.assertEqual(r.returncode,2)

    def test_reject_missing_mpfr_nu(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'trace.txt'
            p.write_text('0 1 0 1 0\n')
            with self.assertRaises(ValueError):load_mpfr_trace(p)

    def test_current_full_independent_result(self):
        # Only data-independent syntax checks unless current P237 source supplied.
        import json
        p=ROOT/'P242_AUDIT_CERT.json'
        x=json.loads(p.read_text())
        self.assertEqual(x['total_independent_positive_pivots'],754)
        self.assertEqual(x['sectors']['even']['independent_mpfr_pivot_count'],377)
        self.assertEqual(x['sectors']['odd']['independent_mpfr_pivot_count'],377)
        self.assertFalse(x['RH_proved'])

if __name__=='__main__':unittest.main()
