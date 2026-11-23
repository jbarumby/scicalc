"""Regression tests for the P92-O variable-width exact rational certificate."""
import unittest
from fractions import Fraction as F
from .certificate import (
    SCHEMA, T, certify, clusters, comparison, pairwise_gate,
    root, examples, hs_square, ldlt_positive, full_delta
)

class P92OTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_scientific_scope(self):
        c=self.data
        self.assertEqual(c["schema"],SCHEMA)
        self.assertEqual(c["case_count"],9)
        self.assertEqual(c["exact_comparison_LDL"],9)
        self.assertEqual(c["full_delta_LDL"],1)
        self.assertGreaterEqual(c["HS_checks"],450)
        for key in ("all_configurations","Preissmann_general","Riemann_hypothesis"):
            self.assertFalse(c[key])

    def test_rational_root(self):
        self.assertEqual(root(F(1,1024)), F(1,32))
        with self.assertRaises(ValueError):
            root(F(2,3))

    def test_block_geometry(self):
        for _,(specs,gaps) in examples().items():
            gs=clusters(specs,gaps)
            for i,g in enumerate(gs):
                self.assertEqual(g["end"]-g["start"],g["w"])
                if i: self.assertEqual(g["start"]-gs[i-1]["end"],gaps[i-1])

    def test_pairwise_three_quarters(self):
        c={r["name"]:r for r in self.data["cases"]}
        self.assertFalse(c["short_spectral_only"]["pairwise_3over4"])
        self.assertFalse(c["anchor_plus_12_microgroups"]["pairwise_3over4"])
        for name in ("six_irregular_not_N","variable_K3","variable_K25","thirty_giant_clusters"):
            self.assertTrue(c[name]["pairwise_3over4"])

    def test_variable_width_outside_prior_N(self):
        r=next(x for x in self.data["cases"] if x["name"]=="six_irregular_not_N")
        self.assertFalse(r["old_P92N"])
        self.assertTrue(r["pairwise_3over4"])

    def test_row_sum_failure_but_ldl_succeeds(self):
        r=next(x for x in self.data["cases"] if x["name"]=="anchor_plus_12_microgroups")
        self.assertGreater(F(r["max_row"]),T)
        self.assertGreater(F(r["ldlt_pivot_lower"]),3)

    def test_hs_pairwise_independently(self):
        s,g=examples()["six_irregular_not_N"]
        gs=clusters(s,g); A=comparison(gs)
        for i in range(len(gs)):
            for j in range(i+1,len(gs)):
                self.assertLessEqual(hs_square(gs,i,j),A[i][j]**2)

    def test_independent_full_delta_positive(self):
        s,g=examples()["short_spectral_only"]
        gs=clusters(s,g); D=full_delta(gs)
        self.assertEqual(len(D),13)
        p=ldlt_positive([[T*(i==j)-D[i][j] for j in range(13)] for i in range(13)])
        self.assertGreater(p,0)

    def test_exact_target_comparison(self):
        self.assertLess(F(5477,5000)**2,F(6,5))
        self.assertGreater(F(157,50)**2/F(3)*F(5477,5000),T)

    def test_rejections(self):
        with self.assertRaises(ValueError): clusters([(1,1)],[])
        with self.assertRaises(ValueError): clusters([(2,1),(2,4)],[1])
        with self.assertRaises(AssertionError): ldlt_positive([[0]])
        with self.assertRaises(ValueError): full_delta(clusters([(100,1)],[]))

if __name__=="__main__":
    unittest.main()
