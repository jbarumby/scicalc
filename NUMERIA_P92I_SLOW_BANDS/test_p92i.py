"""P92I exact regression tests. The all-real theorems are proved in README."""
import unittest
from fractions import Fraction as F
from .certificate import (SCHEMA, ZETA2_UP, RESIDUAL_LOWER,
                          RATIO_ROOT_EXTRA, coefficient, exact_root,
                          root_upper, sample_geometry, delta_matrix,
                          integral_bound, exact_ldlt, verify, certify)

class P92ITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_census_and_scientific_scope(self):
        p=self.data
        self.assertEqual(p["schema"],SCHEMA)
        self.assertEqual((p["rational_regressions"],p["two_sided_LDL_cases"]),(6,4))
        self.assertTrue(p["all_real_bounded_spread_proved_analytically"])
        self.assertTrue(p["irregular_block_comparator_proved"])
        for k in ("arbitrary_configurations_global_residual_proved",
                  "general_Preissmann_Lemma6_proved","Riemann_Hypothesis_proved"):
            self.assertFalse(p[k],k)

    def test_exact_pi_directions_and_factor_four(self):
        self.assertLess(F(242,147),ZETA2_UP)
        self.assertEqual(RESIDUAL_LOWER,F(2725,768))
        self.assertLess(F(109,100)**2,F(6,5))
        self.assertLess(2*ZETA2_UP,RESIDUAL_LOWER)
        self.assertEqual(coefficient(2),2)

    def test_refined_band_strict_residual(self):
        s=RATIO_ROOT_EXTRA
        c=coefficient(s)
        self.assertEqual(F(self.data["refined_ratio"]),s*s)
        self.assertEqual(F(self.data["refined_polynomial_margin"]),F(24,5)-c*c)
        self.assertGreater(F(24,5)-c*c,0)

    def test_local_envelope_entire_fraction_grid(self):
        for a in range(10,21):
            for b in range(10,21):
                aa,bb=F(a,10),F(b,10)
                self.assertLessEqual(aa*bb*(aa-bb)**2,2)

    def test_nearest_neighbor_roots_rational(self):
        x,d,roots,blocks=sample_geometry(
            [[F(1),F(6,5),F(7,5)],[F(10),F(11),F(12)]],(F(500),))
        self.assertEqual(len(blocks),2)
        self.assertEqual([r*r for r in roots],d)
        self.assertEqual(d,[min(abs(x[i]-x[j]) for j in range(len(x)) if j!=i)
                            for i in range(len(x))])

    def test_hs_integral_ceiling_and_rational_root_upper(self):
        for row in self.data["cases"]:
            for pair in row["cross_block_data"]:
                actual=F(pair["actual_HS_squared"])
                ceiling=F(pair["integral_HS_squared_upper"])
                root=F(pair["block_norm_upper"])
                self.assertLessEqual(actual,ceiling)
                self.assertLessEqual(ceiling,root*root)
                self.assertGreaterEqual(root,0)

    def test_all_six_exact_block_gates(self):
        for row in self.data["cases"]:
            self.assertTrue(row["block_test_passes"])
            self.assertLess(F(row["rational_block_norm_upper"]),RESIDUAL_LOWER)
        self.assertEqual(F(self.data["cases"][2]["global_gap_spread"]),676)
        self.assertEqual(F(self.data["cases"][5]["global_gap_spread"]),4225)

    def test_independent_two_sided_ldlt_exact(self):
        self.assertEqual(sum(x["two_sided_exact_ldlt_min_pivots"] is not None
                             for x in self.data["cases"]),4)
        for row in self.data["cases"]:
            p=row["two_sided_exact_ldlt_min_pivots"]
            if p is not None:
                self.assertEqual(len(p),2)
                self.assertTrue(all(F(t)>0 for t in p))

    def test_invalid_input_and_rational_square_root(self):
        self.assertEqual(exact_root(F(49,25)),F(7,5))
        with self.assertRaises(ArithmeticError): exact_root(F(5,4))
        with self.assertRaises(ValueError): coefficient(F(7,2))
        with self.assertRaises(ValueError): root_upper(-1)
        self.assertGreaterEqual(root_upper(F(5,3))**2,F(5,3))
        with self.assertRaises(ValueError): sample_geometry([],())

if __name__=="__main__":
    unittest.main()
