"""P92-H exact rational tests; proof for all real data is in README.md."""
import unittest
from fractions import Fraction as Q
from .certificate import (
    SCHEMA, THRESHOLD_GAP_RATIO, THRESHOLD_SQRT_RATIO, T, TARGET,
    PHI_MAX, ROW_BOUND, phi, point_data, delta_matrix, index_envelope,
    exact_ldlt_min_pivot, verify_example, certify,
)


class P92HLacunaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=certify()

    def test_exact_constants(self):
        self.assertEqual(THRESHOLD_GAP_RATIO,Q(64,25))
        self.assertEqual(THRESHOLD_SQRT_RATIO,Q(8,5))
        self.assertEqual(ROW_BOUND,Q(12793,13824))
        self.assertEqual(TARGET-ROW_BOUND,Q(167,13824))
        self.assertEqual(ROW_BOUND,Q(4,9)+Q(625,1536)+Q(2,27))

    def test_reproduction_and_scope(self):
        cert=self.cert
        self.assertEqual(cert['schema'],SCHEMA)
        self.assertEqual(cert['exact_rational_regression_count'],71)
        self.assertEqual(cert['two_sided_exact_ldlt_cases'],35)
        self.assertTrue(cert['all_lacunary_variable_ratios_proved_analytically'])
        for flag in ('geometric_exact_ratio_required',
                     'all_gap_ratios_above_one_proved_to_meet_15_over_16',
                     'arbitrary_point_configurations_proved',
                     'full_Preissmann_Lemma6_proved','Riemann_Hypothesis_proved'):
            self.assertFalse(cert[flag],flag)

    def test_nearest_neighbor_identity_for_variable_ratios(self):
        for ratios in ((Q(8,5),), (Q(8,5),Q(5,3),Q(3)),
                       (Q(2),Q(10),Q(8,5),Q(9,4),Q(7,4))):
            xs,rr=point_data(ratios)
            self.assertEqual(len(xs),len(rr))
            self.assertEqual(xs[0],0)
            self.assertEqual(rr[0],rr[1])
            for i in range(1,len(rr)-1):
                self.assertGreaterEqual(rr[i+1]/rr[i],Q(8,5))

    def test_scalar_phi_envelope_exact(self):
        self.assertEqual(phi(Q(1,3)),PHI_MAX)
        for d in range(1,18):
            for n in range(d+1):
                s=Q(n,d)
                self.assertLessEqual(phi(s),PHI_MAX)
                self.assertLessEqual(phi(s),s)
        for k in range(1,25):
            self.assertLessEqual(index_envelope(k),PHI_MAX)
            self.assertLessEqual(index_envelope(k),T**k)

    def test_anchor_dominance_and_block_kernel(self):
        for ratios in ((Q(8,5),)*6, (Q(8,5),Q(3),Q(2),Q(9,4),Q(5),Q(8,5))):
            xs,rr=point_data(ratios)
            mat=delta_matrix(xs,rr)
            self.assertEqual(mat[0][1],0)
            for j in range(2,len(xs)):
                self.assertLessEqual(mat[0][j],mat[1][j])
            for i in range(1,len(xs)):
                for j in range(i+1,len(xs)):
                    self.assertLessEqual(mat[i][j],index_envelope(j-i)/2)

    def test_exact_schur_bound_all_regressions(self):
        for case in self.cert['cases']:
            self.assertGreaterEqual(Q(case['margin_to_rational_bound']),0)
            self.assertLessEqual(Q(case['max_exact_row_sum']),ROW_BOUND)
            self.assertEqual(
                ROW_BOUND-Q(case['max_exact_row_sum']),
                Q(case['margin_to_rational_bound']))

    def test_two_sided_ldlt_independent(self):
        for ratios in ((), (Q(8,5),), (Q(8,5),Q(8,5),Q(9,4)),
                       (Q(5),Q(8,5),Q(2),Q(7,4))):
            xs,rr=point_data(ratios)
            mat=delta_matrix(xs,rr)
            self.assertGreater(exact_ldlt_min_pivot(mat,+1),0)
            self.assertGreater(exact_ldlt_min_pivot(mat,-1),0)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): point_data((Q(3,2),))
        with self.assertRaises(ValueError): phi(Q(2))
        with self.assertRaises(ValueError): index_envelope(0)
        with self.assertRaises(ValueError): delta_matrix((Q(0),Q(0)),(Q(1),Q(1)))


if __name__=='__main__':
    unittest.main()
