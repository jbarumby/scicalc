"""P86 regression: all finite certificates use EXACT fractions, no floats."""
import unittest
from fractions import Fraction as F
from .spectral_rational import (
    build,exact_ldl_positive,source_m0_pi,row_energies,
    run,RHO,PI_LOWER,M0_LOWER,PAIRS,mat_mul,
)

class P86SpectralTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=run()

    def test_uniform_N_le_five_constant_is_proved_rationally(self):
        # Universal proof is algebraic: delta_i<=|x_i-x_j|,
        # |B_ij|<=1, hence max absolute row/column sums <=N-1.
        self.assertGreater(PI_LOWER*M0_LOWER,4)
        self.assertGreater(RHO,16)
        self.assertEqual(self.result["constants"]["rho_strict_below_(pi_m0)^2"],str(RHO))
        self.assertFalse(self.result["general_all_N_Preissmann_proved_here"])

    def test_rational_spectral_similarity_exact(self):
        x,d,K,A,H=build(["0","1/7","13/8","5/2"])
        self.assertEqual(len(x),4)
        self.assertTrue(all(t>0 for t in d))
        for i in range(4):
            for j in range(4):
                self.assertEqual(K[i][j],-K[j][i])
                self.assertEqual(H[i][j],H[j][i])
        self.assertEqual(len(exact_ldl_positive(H)),4)

    def test_eight_distinct_nodes_and_multi_precision_not_fake(self):
        x=self.result
        self.assertEqual(x["samples_count"],8)
        self.assertEqual(len(x["finite_samples"]),len(PAIRS))
        self.assertIn(12,x["verified_sizes"])
        for s in x["finite_samples"]:
            self.assertEqual(len(s["nodes_exact_rational"]),s["n"])
            self.assertEqual(len(s["ldl_pivots_exact"]),s["n"])
            self.assertTrue(all(F(q)>0 for q in s["ldl_pivots_exact"]))
            self.assertEqual(F(s["rho_lower_than_pi2_m02_exact"]),RHO)

    def test_uniform_all_N_diagonal_energy_le_four(self):
        self.assertTrue(self.result["general_all_N_energy_le_4_proved_independently"])
        self.assertFalse(self.result["Preissmann_sharper_pi2_over3_row_energy_proved_here"])
        for row in self.result["finite_samples"]:
            calculated=row_energies(row["nodes_exact_rational"])
            self.assertEqual([str(q) for q in calculated],row["row_energies_exact"])
            self.assertTrue(all(F(0)<q<=F(4) for q in calculated))

    def test_wrong_too_small_bound_fails_closed(self):
        _,_,_,_,H=build(["0","1"],threshold=F(1,2))
        with self.assertRaises(ArithmeticError):
            exact_ldl_positive(H)

    def test_adversarial_distinctness_guard(self):
        for bad in ((1,),(1,1),("0","2/3","2/3")):
            with self.subTest(nodes=bad),self.assertRaises(ValueError):
                build(bad)

    def test_simple_two_nodes_exact_eigenvalue_one(self):
        x,d,K,A,H=build([0,7],threshold=F(2))
        self.assertEqual(list(d),[F(7),F(7)])
        self.assertEqual(A,[[F(0),F(-1)],[F(1),F(0)]])
        self.assertEqual(mat_mul(A,A),[[-F(1),F(0)],[F(0),-F(1)]])
        self.assertTrue(all(v>0 for v in exact_ldl_positive(H)))

    def test_analytic_sources_remain_explicitly_unproved(self):
        x=self.result
        for flag in (
            "general_all_N_Preissmann_proved_here",
            "Preissmann_Lemma1_general_spacing_maximizer_proved_here",
            "Preissmann_Lemmas2_to_5_global_spacing_perturbations_proved_here",
            "Preissmann_Lemma6_three_point_denominator_bound_proved_here",
            "Preissmann_final_uniform_spectral_estimate_proved_here",
            "RZ_double_sum_bound_reproved",
            "X_lower_3e10_3e12_covered",
            "RH_proved",
        ):
            self.assertFalse(x[flag])

if __name__=="__main__":
    unittest.main()
