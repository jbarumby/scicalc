"""Independent tests for the substitution and explicit integral tail.

These double-precision examples are diagnostic. The general inequalities
are proved in README.md and must be evaluated in Arb for the actual bound.
"""
import math
import unittest
from fractions import Fraction
from NUMERIA_P76_COEFFICIENTS.table import TABLE

class AlgebraicChecks(unittest.TestCase):
    def test_table(self):
        self.assertEqual(len(TABLE),16)
        self.assertEqual([r.sigma_start for r in TABLE],list(range(16,32)))
        self.assertTrue(all(Fraction(r.log10_a0)>10 for r in TABLE))

    def test_substitution(self):
        for x in (0.1,0.2,0.7,1.0,1.9):
            self.assertAlmostEqual((x**6)**(1/6)*6*x**5,6*x**6,places=12)
        self.assertEqual(2**6,64)

    def test_tail_polynomial(self):
        from scipy.integrate import quad
        for u in (3,10,20):
            value=quad(lambda t:math.exp(-t)*t**3,u,100)[0]
            closed=math.exp(-u)*(u**3+3*u*u+6*u+6)
            self.assertAlmostEqual(value,closed,delta=closed*3e-9)

    def test_tail_sample(self):
        for kap in (.07,.20,.34):
            c=3/(4*kap);C1=.65;C2=.15;L=math.log(3e12)
            for u in (64,80,128,256):
                lhs=u**(1/6)*(1+math.log1p(c*u)/(2*L))**(4/3)
                lhs*= (C1*math.sqrt(c*u)+C2)**(2/3)
                rhs=(1+c/(2*L))**2*(1+C1*math.sqrt(c)+C2)*u**3
                self.assertLessEqual(lhs,rhs)

if __name__=="__main__":unittest.main()
