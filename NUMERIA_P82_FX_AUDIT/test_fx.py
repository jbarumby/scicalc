"""Independent P82 regression tests.

Rational arithmetic proves finite coefficient identities and exact algebra.
The KLN 3.14 mean-value inequality and Ramaré short-range b3=0.605
remain cited theorems; numeric moment test is diagnostic only.
"""
import unittest
from fractions import Fraction as Q
from .fx_audit import (
    mu,lam,divisor_raw,truncated_diagonal,finite_splitting,
    rational_mv_weight_identity,main_audit,
)

class FXAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc=main_audit()

    def test_mobius_basic_values(self):
        known={1:1,2:-1,3:-1,4:0,5:-1,6:1,7:-1,8:0,
               9:0,10:1,12:0,15:1,30:-1,36:0}
        for n,v in known.items():
            with self.subTest(n=n):
                self.assertEqual(mu(n),v)

    def test_coefficients_vanish_below_x(self):
        for X in range(2,18):
            for n in range(1,X+1):
                self.assertEqual(lam(n,X),0)
                if n==1:self.assertEqual(divisor_raw(n,X),1)

    def test_minus_one_is_essential(self):
        # A_X(n) differs from lambda_X(n) at n=1 ONLY.
        for X in (2,5,11):
            n=3*X+1
            a=truncated_diagonal(X,n,3)
            self.assertEqual(Q(a["difference_exact"]),Q(1))
            self.assertEqual(Q(a["raw_diagonal"]),Q(a["lambda_square_sum"])+1)
            self.assertEqual(Q(a["double_sum_truncation"]),Q(a["raw_diagonal"]))

    def test_all_ten_exact_diagonal_cases(self):
        self.assertEqual(len(self.doc["finite_rational_mobius_examples"]),10)
        for v in self.doc["finite_rational_mobius_examples"]:
            self.assertTrue(v["all_equalities_exact"])
            self.assertEqual(Q(v["difference_exact"]),1)

    def test_KLN_m0_weight_algebra(self):
        cases=self.doc["finite_rational_mean_weight_examples"]
        self.assertEqual(len(cases),6)
        for v in cases:
            self.assertEqual(Q(v["direct"]),Q(v["reassembled"]))
            self.assertEqual(Q(v["direct"]),Q(v["term_first_sum"])+Q(v["term_second_sum"]))

    def test_five_x_split_exact(self):
        cases=self.doc["finite_rational_tail_examples"]
        self.assertEqual(len(cases),6)
        self.assertTrue(all(x["combined_satisfied"] for x in cases))

    def test_adversarial_inputs_fail(self):
        for call in (
            lambda:mu(0),lambda:mu(-1),
            lambda:truncated_diagonal(2,1,2),
            lambda:truncated_diagonal(2,8,1),
            lambda:finite_splitting(3,20,0),
            lambda:rational_mv_weight_identity(2,20,2,1,1),
            lambda:rational_mv_weight_identity(2,20,3,-1,1),
        ):
            with self.assertRaises(ValueError):
                call()

    def test_external_gates_not_misrepresented(self):
        d=self.doc
        self.assertIn("pi*m0",d["analytic_integration"]["inequality"])
        self.assertIn("5X",d["tail_splitting"]["long_range"])
        self.assertIn("0.605",d["tail_splitting"]["b3_short_range_imported"])
        self.assertTrue(d["status"]["CS_FX_inequality_derived_conditionally"])
        self.assertFalse(d["status"]["KLN_mean_value_theorem_independently_proved"])
        self.assertFalse(d["status"]["KLN_b3_short_range_independently_proved"])
        self.assertFalse(d["status"]["RZ_S_epsilon_theorems_independently_proved"])
        self.assertFalse(d["status"]["Riemann_hypothesis_proved"])
        self.assertEqual(d["final_constants"]["for_CS_corollary_1"],
                         "X=h*T, h=1, T>=3e12 (CS proof of Corollary 1)")

    def test_independent_finite_moment_diagnostic(self):
        # Exact formula for finite oscillatory integral, evaluated numerically.
        # It is a sanity check on signs/normalization and NOT a theorem proof.
        import mpmath as mp
        mp.mp.dps=65
        m0=mp.sqrt(1+mp.mpf(2)/3*mp.sqrt(mp.mpf(6)/5))
        for X,N,T in ((2,35,mp.mpf('0.7')),
                      (4,42,mp.mpf('3.0')),
                      (7,48,mp.mpf('10'))):
            sigma=mp.mpf('1.5')
            u={n:mp.mpf(lam(n,X))/mp.power(n,sigma)
               for n in range(1,N+1) if lam(n,X)!=0}
            diagonal=T*mp.fsum(t*t for t in u.values())
            cross=mp.mpf(0)
            for m in u:
                for n in u:
                    if m<n:
                        d=mp.log(mp.mpf(n)/m)
                        cross+=2*u[m]*u[n]*mp.sin(T*d)/d
            integral=diagonal+cross
            KLN_upper=mp.fsum(u[n]**2*(T+mp.pi*m0*(n+1)) for n in u)
            with self.subTest(X=X,T=str(T)):
                self.assertTrue(0<=integral<KLN_upper)

if __name__=="__main__":
    unittest.main()
