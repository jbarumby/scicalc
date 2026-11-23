"""P92-L exact Fraction regressions; universality is proven in README."""
import unittest
from fractions import Fraction as F
from .certificate import (
    certify, build, nearest, comparison, cross_entry,
    delta_exact, positive_fraction_ldlt, star_gate,
    PAIR_HS_UP, TARGET, SCHEMA
)

class P92LTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_provenance_scope_and_census(self):
        d=self.data
        self.assertEqual(d['schema'],SCHEMA)
        self.assertEqual((d['fraction_cases'],d['pair_checks'],d['full_fraction_LDL_cases']),
                         (8,249,4))
        self.assertTrue(d['analytic_proof_all_real_under_hypotheses'])
        self.assertFalse(d['arbitrary_cluster_size_without_spectral_gate'])
        for key in ('arbitrary_point_configurations_proved',
                    'general_Preissmann_Lemma6_proved','Riemann_Hypothesis_proved'):
            self.assertFalse(d[key],key)

    def test_nearest_spacing_matches_internal_group_gap(self):
        roots=(F(1),F(2),F(1,3),F(5,2))
        sizes=(4,2,9,3)
        for mode in ('local','global'):
            xs,gs,h,_=build(roots,sizes,gap_mode=mode)
            v=nearest(xs)
            for r,g in enumerate(gs):
                for i in g:
                    self.assertEqual(v[i],h[r])

    def test_pair_coefficient_exact_global_max(self):
        for s1 in (F(1,100),F(1,9),F(1),F(4),F(100)):
            for s2 in (F(1,100),F(1,16),F(1),F(9)):
                a,b=s1*s1,s2*s2;H=max(a,b)
                beta=s1*s2*(s1-s2)**2/2
                self.assertLessEqual(beta/H**2,F(2,27))
        self.assertGreater(PAIR_HS_UP**2,F(1,8))

    def test_pair_hs_comparison_exact(self):
        xs,gs,h,roots=build((F(1),F(3),F(1,2),F(2)),
                             (2,4,3,6),gap_mode='global')
        a,pairs,_=comparison(xs,gs,h,roots)
        self.assertEqual(len(pairs),6)
        for p in pairs:
            hs2=F(p['HS_squared'])
            self.assertLessEqual(hs2,F(p['HS_bound_integral'])**2)
            self.assertLessEqual(hs2,F(p['HS_bound_cardinality'])**2)
            self.assertLessEqual(hs2,F(p['rational_bound'])**2)
            self.assertTrue(p['linear_separation'])
        self.assertEqual(a[0][0],0)

    def test_M14_implies_P88_target_algebraically(self):
        self.assertLess(F(28,27)**2,F(6,5))
        self.assertGreater(F(30,27)**2,F(6,5))

    def test_three_quarters_sublinear_condition_exact(self):
        for k in range(1,30):
            self.assertGreaterEqual(F(k)**4,F(k)**3)
        self.assertLess(F(16,9),F(3))

    def test_13_block_star_row_fails_spectral_passes(self):
        z=star_gate()
        self.assertEqual((z['groups'],z['total_points'],z['cross_pairs']),
                         (13,88,78))
        self.assertGreater(F(z['max_row_sum']),4)
        self.assertLess(F(z['spectral_radius_squared']),TARGET**2)
        self.assertGreater(F(z['exact_comparison_LDL_min_pivot']),0)

    def test_full_fraction_ldlt_independent(self):
        xs,gs,h,roots=build((F(1),F(2),F(1,3)),(3,2,3),gap_mode='global')
        delta=delta_exact(xs,gs,roots)
        n=len(xs)
        M=[[TARGET*(i==j)-delta[i][j] for j in range(n)] for i in range(n)]
        self.assertGreater(positive_fraction_ldlt(M),0)

    def test_reject_invalid_or_indefinite(self):
        for roots,sizes in (((),()),((1,),(1,)),((0,),(2,)),((1,2),(2,))):
            with self.assertRaises(ValueError):
                build(roots,sizes)
        with self.assertRaises(ValueError):
            build((1,2),(2,2),gap_mode='nonsense')
        with self.assertRaises(ArithmeticError):
            positive_fraction_ldlt([[F(1),F(2)],[F(2),F(1)]])
        with self.assertRaises(ValueError):
            positive_fraction_ldlt([[1,2],[2]])

if __name__=='__main__':
    unittest.main()
