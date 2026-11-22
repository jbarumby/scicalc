"""P81 regression tests. Theorems by Ramaré are imported, not proved."""
import unittest
from fractions import Fraction as F
from .transfer import audit,check_positive_bound,SCHEMA

class TransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=audit()

    def test_full_16x3x2_coverage(self):
        r=self.result
        self.assertEqual(r["schema"],SCHEMA)
        self.assertEqual(r["total_branch_zone_checks"],96)
        self.assertEqual([z["row"] for z in r["rows"]],list(range(1,17)))
        self.assertTrue(r["all_checks_pass"])
        self.assertGreater(F(r["smallest_rational_b4_slack"]),0)

    def test_each_branch_uses_same_authorized_public_bound(self):
        for row in self.result["rows"]:
            for zone in row["zones"]:
                self.assertEqual(len(zone["branches"]),2)
                eps2,eps1=zone["branches"]
                self.assertEqual(eps2["branch"],"epsilon_2")
                self.assertEqual(eps1["branch"],"epsilon_1")
                self.assertEqual(F(eps2["epsilon_zone_max"]),2*F(eps1["epsilon_zone_max"]))
                self.assertEqual(F(eps2["source_S_cap"]),F(eps1["source_S_cap"]))
                self.assertEqual(F(eps2["source_epsilon_cap"]),F(eps1["source_epsilon_cap"]))
                self.assertLessEqual(F(eps2["epsilon_zone_max"]),F(eps2["source_epsilon_cap"]))
                self.assertGreater(F(eps1["slack_under_b4"]),0)
                self.assertGreater(F(eps2["slack_under_b4"]),0)
                self.assertLess(F(eps2["slack_under_b4"]),F(eps1["slack_under_b4"]))
                self.assertGreater(F(row["b4_from_P80"]),F(eps2["zeta_product_bound"]))

    def test_expected_v3_improvement_preserved(self):
        r=self.result["rows"][9]
        self.assertEqual(F(r["b4_from_P80"]),F(139,250))
        self.assertEqual(F(r["b4_previous_P78"]),F(557,1000))
        self.assertGreater(F(r["b4_from_P80"]),F(r["zones"][0]["branches"][0]["zeta_product_bound"]))

    def test_counterexample_to_bad_b4_is_rejected(self):
        with self.assertRaises(ArithmeticError):
            check_positive_bound(F(1,50),F(5501,10000),F(1,2))
        with self.assertRaises(ArithmeticError):
            check_positive_bound(F(1,2),F(1,2),F(1))
        with self.assertRaises(ArithmeticError):
            check_positive_bound(F(1,100),F(1,2),F(0))

    def test_no_global_claim(self):
        r=self.result
        self.assertFalse(r["external_S_epsilon_theorems_proved_here"])
        self.assertFalse(r["original_CS_Lemma6_full_X_range_proved_here"])
        self.assertFalse(r["independent_KLN_F_X_lemma_proved_here"])
        self.assertFalse(r["global_RH_proved"])
        self.assertFalse(r["Suzuki_new_window_proved_here"])

    def test_C4_algebra_has_5_power(self):
        c=self.result["checked_analytic_implications"]["eps1_C4"]
        self.assertIn("5^(-eps1)",c)
        self.assertIn("e^(-delta)",c)
        self.assertIn("logX",c)

if __name__=="__main__":
    unittest.main()
