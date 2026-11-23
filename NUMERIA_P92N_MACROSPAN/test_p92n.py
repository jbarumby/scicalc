"""P92-N regression suite in exact rationals; analytic proof in README."""
import unittest
from fractions import Fraction as F
from .certificate import (
    C,Q,T,SCHEMA,certify,build,full_delta,
    exact_ldl,rational_threshold,large_macro,audit
)

class P92NTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_scope_and_census(self):
        d=self.data
        self.assertEqual(d["schema"],SCHEMA)
        self.assertEqual(d["rational_regressions"],6)
        self.assertEqual(d["exact_pair_hs_checks"],56)
        self.assertEqual(d["full_exact_ldl_count"],3)
        self.assertTrue(d["analytic_theorem_for_all_real_macroscale_configurations"])
        for k in ("unrestricted_point_configurations_proved",
                  "full_Preissmann_Lemma6_proved","Riemann_Hypothesis_proved"):
            self.assertFalse(d[k],k)

    def test_machin_fraction_threshold(self):
        self.assertEqual(rational_threshold(),F(619520,173901))
        self.assertLess(rational_threshold(),F(18,5))
        self.assertGreater(Q**2,F(1,2))
        self.assertGreater(F(294)*13**4,F(125)*16**4)

    def test_exact_nearest_neighbor_geometry(self):
        g,x=build([F(1,2),F(1,3),F(1,5),F(1)],
                  [5,10,26,2],materialize=True)
        self.assertEqual([min(abs(x[i]-x[j]) for j in range(len(x)) if i!=j)
                          for i in range(len(x))],
                         [group[1] for group in g for _ in range(group[2])])
        self.assertTrue(all(group[4]-group[3]<=1 for group in g))

    def test_exact_pair_bounds_for_all_regressions(self):
        for case in self.data["cases"]:
            self.assertGreater(F(case["comparison_ldl_min_pivot"]),0)
            self.assertLess(F(case["max_comparison_row_sum"]),T)
            self.assertGreaterEqual(case["pair_hs_checks"],1)

    def test_full_fraction_ldl(self):
        for case in self.data["cases"]:
            if case["full_ldl_min_pivot"] is not None:
                self.assertGreater(F(case["full_ldl_min_pivot"]),0)
        g,_=build([F(1,2),F(1,3)],[3,5])
        D=full_delta(g)
        self.assertGreater(exact_ldl([[T*(i==j)-D[i][j] for j in range(8)] for i in range(8)]),0)

    def test_same_step_means_zero_delta(self):
        g,_=build([F(1,2)]*4,[5]*4)
        M=full_delta(g)
        self.assertTrue(all(x==0 for row in M for x in row))

    def test_40_groups_22180_points_without_full_matrix(self):
        case=large_macro()
        self.assertEqual(case["groups"],40)
        self.assertEqual(case["points"],22180)
        self.assertEqual(case["distinct_scales"],40)
        self.assertEqual(F(case["end_to_end_gap"]),F(559,8))
        self.assertTrue(case["P92M_separation_violated"])
        self.assertFalse(case["full_matrix_constructed"])
        self.assertLess(F(case["macro_comparison_row_sum_exact"]),rational_threshold())

    def test_invalid_physical_spans(self):
        with self.assertRaises(AssertionError):
            build([F(1)], [3])
        with self.assertRaises(AssertionError):
            build([F(1,3)], [1])

if __name__=="__main__":
    unittest.main()
