"""P92K: nine precise Fraction checks for multicluster theorem regression."""
import unittest
from fractions import Fraction as F
from .certificate import (SCHEMA, BLOCK_FACTOR, TARGET, ZERO,
                          certify, construct, nearest, matrix, audit_case)
from NUMERIA_P92J_MULTIANCHOR.certificate import ldlt_positive

class P92KTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=certify()

    def test_01_counts(self):
        d=self.d
        self.assertEqual((d["fraction_configurations"],d["pairwise_interactions"],
                          d["exact_full_LDL"],d["exact_comparison_LDL"]), (11,760,4,11))
        self.assertEqual(len(d["cases"]),11)

    def test_02_scope(self):
        d=self.d
        self.assertEqual(d["schema"],SCHEMA)
        self.assertTrue(d["all_real_parameters_with_quadratic_separation_proved"])
        for key in ("arbitrary_configurations_proved",
                    "general_Preissmann_Lemma6_proved","Riemann_Hypothesis_proved"):
            self.assertFalse(d[key],key)

    def test_03_real_nearest_neighbor_formula(self):
        roots=(F(1,2),F(3),F(1),F(7,5))
        x,g,h,s=construct(roots,(2,3,4,2))
        self.assertEqual(nearest(x),tuple(h[r] for r,group in enumerate(g) for _ in group))

    def test_04_exact_cross_geometry(self):
        roots=(F(1),F(5),F(1,7),F(3))
        x,groups,h,s=construct(roots,(2,3,4,2))
        for r in range(len(groups)):
            for t in range(r+1,len(groups)):
                D=x[groups[t][0]]-x[groups[r][-1]]
                self.assertGreaterEqual(D,(t-r)**2*max(h[r],h[t]))

    def test_05_zero_intrablocks_and_nonnegative_offblocks(self):
        roots=(F(1),F(2),F(1,2))
        x,groups,h,s=construct(roots,(2,3,2))
        roots_node=tuple(s[r] for r,g in enumerate(groups) for _ in g)
        delta=matrix(x,roots_node)
        for g in groups:
            self.assertTrue(all(delta[i][j]==0 for i in g for j in g))
        self.assertTrue(all(delta[i][j]==delta[j][i]>=0 for i in range(len(x)) for j in range(len(x))))

    def test_06_cross_pair_hs_enclosures(self):
        for case in self.d["cases"]:
            self.assertEqual(case["pair_count"],case["group_count"]*(case["group_count"]-1)//2)
            for pair in case["pairs"]:
                self.assertGreaterEqual(F(pair["rho"]),(pair["s"]-pair["r"])**2)
                self.assertLessEqual(F(pair["hs_square"]),F(pair["hs_square_ceiling"]))
                self.assertEqual(F(pair["hs_square_ceiling"]),F(1,8)/F(pair["rho"])**2)

    def test_07_real_polynomial_monotonicity_gate(self):
        for n in (1,3,7,19,100):
            for i in range(n+1):
                u=F(i,n)
                self.assertLessEqual((1-u)**4*(3+8*u*u),3)
                self.assertGreater(12-16*u+48*u*u,0)

    def test_08_rational_schur_comparison(self):
        self.assertGreater(BLOCK_FACTOR**2,F(1,8))
        self.assertTrue(all(F(c["max_rational_comparison_row_sum"])<F(10,7) for c in self.d["cases"]))
        self.assertTrue(ldlt_positive([[TARGET,F(-1,8)],[F(-1,8),TARGET]]))

    def test_09_reject_invalid_clusters_and_run_mixed(self):
        with self.assertRaises(ValueError):
            construct((F(1),),(1,))
        with self.assertRaises(ValueError):
            construct((F(0),),(2,))
        with self.assertRaises(ValueError):
            construct((F(1),F(2)),(2,))
        result=audit_case((F(1),F(5,3),F(1,7)),(3,2,4),True)
        self.assertEqual(result["group_count"],3)
        self.assertTrue(result["full_matrix_fraction_ldlt"])

if __name__=="__main__":
    unittest.main()
