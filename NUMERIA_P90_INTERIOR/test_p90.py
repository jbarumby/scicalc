"""P90 exact Fraction regression: two universal mathematical claims are proved in README."""
import unittest
from fractions import Fraction as F
from .nearby_sharp import (
  zeta2_lower,H,uniform_budget_identity,chord_knot_values,
  one_interior_poly,mid_scaled,D_lower,budget_instance,
  sharp_pair_direct,proof_ledger,SCHEMA
)
from .variational import one_point_profile,compare_profiles

class P90Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=proof_ledger()

    def test_knots_and_exact_zetavalues(self):
        k=chord_knot_values()
        self.assertEqual(k["S_at_zero"],"Z")
        self.assertEqual(k["S_at_half"],"2*Z-11/4")
        self.assertEqual(k["S_at_one"],"2*Z-3")
        self.assertGreater(zeta2_lower(),F(3,2))
        self.assertEqual(H(3,2),F(49,36))

    def test_all_uniform_interior_equalities(self):
        self.assertEqual(len(self.d["uniform_all_m_closed_forms"]),40)
        for j in range(1,41):
            q=uniform_budget_identity(j)
            self.assertEqual(F(q["interior_exact"]),2*F(q["D_exact"]))
            self.assertEqual(F(q["two_sided_budget_exact"]),F(q["interior_exact"]))

    def test_symbolic_terminal_factorization_all_t_in_0_half(self):
        for k in range(1,51):
            w=one_interior_poly(F(k,100))
            self.assertTrue(w["factorization_exact"])
            self.assertGreaterEqual(F(w["positive_polynomial_exact"]),0)
        self.assertEqual(F(one_interior_poly(F(1,2))["positive_polynomial_exact"]),0)

    def test_sharp_all_real_near_pairs_are_distinguished(self):
        self.assertEqual(self.d["schema"],SCHEMA)
        self.assertTrue(self.d["universal_adjacent_pairs_proved_in_README"])
        self.assertTrue(self.d["universal_one_interior_pairs_proved_in_README"])
        self.assertEqual(self.d["exact_nearby_pair_count"],96)
        self.assertEqual(sum(x["count"] for x in self.d["exact_nearby_pair_cases"]),96)
        for sample in self.d["exact_nearby_pair_cases"]:
            for pair in sample["nearby_pairs"]:
                self.assertIn(pair["interior_count"],(0,1))
                self.assertTrue(pair["all_real_theorem_covered"])
                self.assertGreater(F(pair["strict_pi_lower_margin"]),0)
                self.assertEqual(F(pair["sharp_RHS_pi_lower_exact"])-
                                 F(pair["sharp_LHS_exact"]),
                                 F(pair["strict_pi_lower_margin"]))

    def test_nonuniform_multiinterior_sufficient_budget_is_rational(self):
        self.assertEqual(len(self.d["interior_sufficient_budget_finite_checks"]),8)
        for q in self.d["interior_sufficient_budget_finite_checks"]:
            self.assertTrue(q["strict_for_selected_example"])
            self.assertEqual(
                F(q["D_endpoint_lower_sum"])-F(q["interior_exact"]),
                F(q["certified_positive_margin"]))
            self.assertGreater(F(q["certified_positive_margin"]),0)

    def test_variational_derivative_nonmonotonic(self):
        x=compare_profiles()
        self.assertTrue(x["one_interior_derivative_changes_sign"])
        self.assertEqual(x["stationary_positions"],["1/3","2/3"])
        self.assertEqual(x["cusp_maximum_position"],"1/2")
        self.assertFalse(x["full_Preissmann_Lemma5_variational_objective_reconstructed"])
        self.assertEqual(F(one_point_profile(F(1,2))["kernel_exact"]),F(8))
        self.assertEqual(F(one_point_profile(F(1,3))["kernel_exact"]),F(27,4))
        self.assertGreater(F(one_point_profile(F(1,4))["kernel_exact"]),
                           F(one_point_profile(F(1,3))["kernel_exact"]))

    def test_source_integrity_and_scope(self):
        d=self.d
        self.assertEqual(d["uploaded_P87"]["ledger_gates_reported"],24)
        self.assertEqual(d["uploaded_P87"]["ledger_gates_failed_reported"],0)
        self.assertTrue(d["uploaded_P87"]["ledger_B5_wrong_upper_pi_and_float_margins"])
        self.assertTrue(d["uploaded_P87"]["cross_control_imported_in_upload"])
        for flag in (
            "arbitrary_interior_count_Lemma6_proved",
            "full_Preissmann_Lemma5_variations_reconstructed",
            "full_Preissmann_Lemma6_independently_reproved",
            "m0_without_external_Lemma6_proved",
            "RZ_S_epsilon_reproved","RH_proved",
        ):
            self.assertFalse(d[flag])

    def test_input_guards(self):
        for x in (F(0),F(1),F(-1)):
            with self.assertRaises(ValueError):one_point_profile(x)
        for x in (F(0),F(3,4)):
            with self.assertRaises(ValueError):one_interior_poly(x)
        for g in ([],[F(1),F(0)],[F(-1),F(1)]):
            with self.assertRaises(ValueError):mid_scaled(g)
        for s in (F(0),F(2)):
            with self.assertRaises(ValueError):D_lower(s)
        with self.assertRaises(ValueError):sharp_pair_direct([0,1,2],1,0)
        with self.assertRaises(ValueError):uniform_budget_identity(0)

if __name__=="__main__":unittest.main()
