"""FLINT/Arb-only tests. Non-finite integrals fail closed."""
import unittest
from flint import arb,ctx
from NUMERIA_P77_B1.audit_b1 import audit,tail_upper
from NUMERIA_P76_COEFFICIENTS.audit_b23 import constants,q

class ArbChecks(unittest.TestCase):
    def test_elementary_tail_bound(self):
        with ctx.workprec(192):
            a=constants();L=a["logh"];pi=a["pi"];m0=a["m0"]
            c1=6/pi**2+q("1.048")/L
            c2=pi*m0*q(".62")/L+6*m0/(pi*3*10**12)+pi*m0*q("1.048")/(arb(3*10**12)*L)
            for k in (".07",".20",".34"):
                r=tail_upper(3/(4*q(k)),c1,c2,L)
                self.assertTrue(r>0 and r<arb("1e-15"))

    def test_first_row(self):
        result=audit(192,one=1)
        self.assertTrue(result["tested_rows_passed"])
        self.assertFalse(result["full_16_B1_verified"])
        self.assertEqual(result["rows"][0]["index"],1)

    def test_bad_precision(self):
        with self.assertRaises(ValueError):audit(64,1)

if __name__=="__main__":unittest.main()
