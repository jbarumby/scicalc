"""P92-J: independent exact arithmetic checks; the continuum proof is in README."""
import unittest
from fractions import Fraction as F
from itertools import product
from .certificate import (
    SCHEMA, BOUND_SQUARE, EASY_BOUND_SQUARE, P88_RESIDUAL_LOWER,
    audit_case, certify, coordinates, delta_matrix_independent,
    ldlt_positive, rectangular_weights, schur_core_reduce,
    sqrt_upper_envelope,
)

class P92JTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=certify()

    def test_scopes_and_census(self):
        a=self.result
        self.assertEqual(a['schema'],SCHEMA)
        self.assertEqual((a['fraction_regressions'],a['reduced_exact_ldlt_cases'],
                          a['independent_full_two_sided_ldlt_cases']), (43,43,28))
        self.assertTrue(a['multi_anchor_family_all_real_parameters_proved'])
        self.assertTrue(a['general_core_leaf_Schur_reduction_proved'])
        for k in ('arbitrary_point_configurations_proved','all_block_partitions_pass_Schur',
                  'full_Preissmann_Lemma6_proved','Riemann_Hypothesis_proved'):
            self.assertFalse(a[k],k)

    def test_nearest_neighbor_geometry(self):
        for s,(r,m) in product((F(1),F(6,5),F(2),F(7)),((1,2),(2,5),(6,4))):
            x,d=coordinates(s,r,m)
            self.assertEqual(len(x),r+m)
            self.assertEqual(d,(s*s,)*r+(F(1),)*m)
            self.assertEqual([x[j+1]-x[j] for j in range(r-1)], [s*s]*(r-1))
            self.assertEqual([x[j+1]-x[j] for j in range(r,len(x)-1)],[F(1)]*(m-1))

    def test_direct_matrix_equals_bipartite_weights(self):
        for s,r,m in ((F(6,5),3,4),(F(3),4,6),(F(1),5,3)):
            matrix=delta_matrix_independent(s,r,m)
            w=rectangular_weights(s,r,m)
            for i in range(r+m):
                for j in range(r+m):
                    if (i<r and j<r) or (i>=r and j>=r):
                        self.assertEqual(matrix[i][j],0)
                    if i<r<=j:
                        self.assertEqual(matrix[i][j],w[r-1-i][j-r])
                        self.assertEqual(matrix[j][i],w[r-1-i][j-r])

    def test_elementary_global_envelope(self):
        self.assertLess(BOUND_SQUARE,EASY_BOUND_SQUARE)
        for n in (1,2,5,11,50):
            for j in range(n+1):
                u=F(j,n)
                self.assertLessEqual((1-u)**4*(1+3*u*u),1)
                self.assertGreater(18*(u-F(1,6))**2+F(7,2),0)
        for s in (F(1),F(11,10),F(3,2),F(4),F(100)):
            self.assertLessEqual(sqrt_upper_envelope(s),BOUND_SQUARE)

    def test_all_fraction_samples(self):
        self.assertEqual(len(self.result['records']),43)
        for row in self.result['records']:
            self.assertLessEqual(F(row['hs_square_exact']),F(row['analytic_rational_envelope']))
            self.assertLessEqual(F(row['analytic_rational_envelope']),BOUND_SQUARE)
            self.assertTrue(row['exact_reduced_schur_ldlt_positive'])
            self.assertEqual(len(row['reduced_schur_pivot_numerator_bits']),row['anchors'])

    def test_generic_core_leaf_schur_exact(self):
        # Nonzero core: an arbitrary nonnegative symmetric 2x2 core and 3 leaves.
        core=[[F(0),F(1,10)],[F(1,10),F(0)]]
        w=[[F(1,12),F(1,10),F(1,8)], [F(1,20),F(1,14),F(1,6)]]
        t=F(1,2)
        pivots=schur_core_reduce(core,w,t)
        self.assertEqual(len(pivots),2)
        z=F(0)
        full=[[core[i][j] if i<2 and j<2 else (
             w[i][j-2] if i<2 else (w[j][i-2] if j<2 else z))
               for j in range(5)] for i in range(5)]
        full_l=[[t*(i==j)-full[i][j] for j in range(5)] for i in range(5)]
        self.assertEqual(len(ldlt_positive(full_l)),5)

    def test_284_points_row_test_fails_spectral_passes(self):
        stress=self.result['records'][-1]
        self.assertEqual((stress['anchors'],stress['cluster_points']),(4,280))
        self.assertEqual(F(stress['row_sum_lower_by_exact_integral']),F(7875,2144))
        self.assertGreater(F(stress['row_sum_lower_by_exact_integral']),P88_RESIDUAL_LOWER)
        self.assertGreater(F(stress['nearest_anchor_row_sum_exact']),P88_RESIDUAL_LOWER)
        self.assertLess(F(stress['hs_square_exact']),BOUND_SQUARE)
        self.assertTrue(stress['exact_reduced_schur_ldlt_positive'])

    def test_invalid_inputs(self):
        for bad in ((F(1,2),1,2),(F(1),0,2),(F(1),1,1)):
            with self.assertRaises(ValueError): coordinates(*bad)
        with self.assertRaises(ValueError): sqrt_upper_envelope(F(0))
        with self.assertRaises(ArithmeticError): ldlt_positive([[F(0)]])
        with self.assertRaises(ValueError): schur_core_reduce([[F(0)]],[[F(1)]],F(0))


if __name__=='__main__':
    unittest.main()
