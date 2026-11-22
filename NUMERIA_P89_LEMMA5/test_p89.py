"""Exact rational P89 regression tests (not a substitute for Shan theorem proof)."""
import unittest
from fractions import Fraction as F
from .exterior import kernel,majorant,exact_parts,audit,source_errata,SCHEMA

class P89Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=audit()

    def test_kernel_exact_and_guards(self):
        self.assertEqual(kernel(F(1),F(1)),F(1,4))
        self.assertEqual(kernel(F(2),F(3)),F(1,100))
        for t,h in ((0,1),(1,0),(-1,2)):
            with self.assertRaises(ValueError):kernel(t,h)

    def test_monotonic_g_kernel_for_rational_samples(self):
        # The analytic proof works for ALL g,h>0, integers j>=1:
        # g*kernel(jg,h)=1/(g*j^2*(jg+h)^2), decreasing in g.
        for h in (F(1,7),F(2),F(101)):
            for j in (1,2,3,10):
                for g in (F(1,100),F(1),F(31)):
                    self.assertGreater(g*kernel(j*g,h),
                                       (g+F(1,7))*kernel(j*(g+F(1,7)),h))

    def test_exact_rational_tail_majorant(self):
        for d,h in ((F(1),F(2)),(F(1,300),F(1)),(F(7,5),F(1000))):
            exact,tail,upper=majorant(d,h,24)
            self.assertEqual(exact+tail,upper)
            self.assertEqual(tail,F(1,3*(24*d)**3))
            self.assertGreater(tail,0)

    def test_actual_exterior_three_region_decomposition(self):
        self.assertEqual(self.cert["schema"],SCHEMA)
        self.assertGreaterEqual(self.cert["finite_checked_pairs"],25)
        for group in self.cert["finite_sets"]:
            self.assertEqual(group["pairs_checked"],len(group["records"]))
            for r in group["records"]:
                parts=[F(r[x]) for x in
                      ("left_exact","middle_exact","right_exact")]
                self.assertEqual(sum(parts,F(0)),F(r["total_exact"]))
                self.assertLessEqual(F(r["left_exact"]),F(r["left_bound_rational"]))
                self.assertLessEqual(F(r["right_exact"]),F(r["right_bound_rational"]))
                self.assertGreaterEqual(F(r["left_gap_rational"]),0)
                self.assertGreaterEqual(F(r["right_gap_rational"]),0)

    def test_degenerate_inputs(self):
        with self.assertRaises(ValueError):exact_parts([0,0,1],0,1)
        with self.assertRaises(ValueError):exact_parts([0,1,2],0,0)
        with self.assertRaises(ValueError):exact_parts([0,1,2],1,0)
        with self.assertRaises(ValueError):exact_parts([0,1,2],0,3)
        with self.assertRaises(ValueError):majorant(0,1)
        with self.assertRaises(ValueError):majorant(1,0)
        with self.assertRaises(ValueError):majorant(1,1,0)

    def test_laplace_correction_and_pi_logic_gate(self):
        z=source_errata()
        self.assertTrue(z["laplace_original_missing_factor_one_over_t"])
        self.assertTrue(z["laplace_original_reversed_exponential_inequality"])
        self.assertTrue(z["B5_original_used_upper_pi_bound_not_sufficient"])
        self.assertGreater(F(z["B5_correct_use_lower_pi_bound"]),F(157,50))
        # Known exact rational witness Q below true lower pi bound squared.
        Q=F(3)
        self.assertLess(Q,F(z["B5_correct_Q_bound_lower"]))

    def test_scope_and_no_false_lemma5_6_claims(self):
        d=self.cert
        self.assertTrue(d["exterior_bound_proved_for_all_finite_real_nodes"])
        for k in ("interior_middle_sharp_Lemma6_proved",
                  "full_Lemma5_variational_proof_reproduced",
                  "full_Lemma6_preissmann_reproved",
                  "RZ_double_sum_source_extended_below_3e12",
                  "RH_proved"):
            self.assertFalse(d[k])

if __name__=="__main__":unittest.main()
