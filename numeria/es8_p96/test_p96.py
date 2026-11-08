"""NUMERIA ES8-P96 tests: no external P89/P95 artifacts are used."""
import ast
from fractions import Fraction as F
from pathlib import Path
import unittest

SRC = Path(__file__).with_name("W12_arb_direto_P96.py").read_text()

class Structure(unittest.TestCase):
    def test_valid_python(self):
        compile(SRC, "W12_arb_direto_P96.py", "exec")
    def test_nonanalytic_guard(self):
        f=SRC[SRC.index("    def outer_integrand"):SRC.index("    @staticmethod\n    def geometric_panels")]
        self.assertLess(f.index("return acb_nan()"),f.index("self.correlation(t + delta)"))
    def test_real_nonnegative_bump(self):
        self.assertIn("e_inv / 2 + arb(0, (e_inv / 2).upper())",SRC)
        self.assertIn("acb(z.real)",SRC)
    def test_exact_extremes(self):
        self.assertIn("exact_real_bounds",SRC)
        self.assertIn(".fmpq()",SRC)
    def test_small_t_error_fraction(self):
        eps=F(1,65536); h=F(2,5)
        self.assertEqual(eps*eps/h,F(5,8589934592))
    def test_support_trimming_rational(self):
        h=F(2,5)
        for a in range(-16,17):
            for b in range(a,17):
                lo=max(-h,-h-F(b,20));hi=min(h,h-F(a,20))
                for si in range(a,b+1):
                    for yi in range(-8,9):
                        s=F(si,20); y=F(yi,20)
                        if -h<y<h and -h<y+s<h:
                            self.assertTrue(lo<=y<=hi)
    def test_no_prior_data_read(self):
        self.assertNotIn("RESULTADO_ES8P89.json",SRC)
        self.assertNotIn("MATRIZ_WEIL_8x8_INTERVALOS_CERTIFICADOS.json",SRC)

if __name__=="__main__":
    unittest.main(verbosity=2)
