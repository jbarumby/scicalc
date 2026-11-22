"""Exact rational multicluster audit, no numerical eigenvalue dependency."""
import unittest
from fractions import Fraction as F
from .certificate import (
    SCHEMA, certify, block_certificate, build_uniform_clusters,
    two_cluster_coefficient, pair_hs_square_exact,
    pair_hs_integral_ceiling_square, sqrt_rational_upper,
)


def exact_delta_matrix(ss, sizes, multipliers):
    """Independent construction from actual nodes and nearest-neighbour distances."""
    ss,sizes,hs,groups=build_uniform_clusters(ss,sizes,multipliers)
    positions=[x for group in groups for x in group]
    cluster_index=[r for r,group in enumerate(groups) for x in group]
    n=len(positions)
    gaps=[]
    for i,x in enumerate(positions):
        gaps.append(min(abs(x-y) for j,y in enumerate(positions) if j!=i))
    selfverify=[hs[cluster_index[i]]==gaps[i] for i in range(n)]
    if not all(selfverify):
        raise ArithmeticError('nearest-neighbour gap mismatch')
    a=[[F(0) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            d1,d2=ss[cluster_index[i]],ss[cluster_index[j]]
            val=d1*d2*(d1-d2)**2/(2*(positions[i]-positions[j])**2)
            a[i][j]=a[j][i]=val
    return a,cluster_index


def exact_spd_ldlt(a):
    """Exact LDL^T positivity (not needed in proof), independent regression."""
    n=len(a)
    L=[[F(0)]*n for _ in range(n)]
    d=[]
    for i in range(n):
        L[i][i]=F(1)
        cur=a[i][i]-sum((L[i][k]**2*d[k] for k in range(i)),F(0))
        if cur<=0:
            return False
        d.append(cur)
        for j in range(i+1,n):
            L[j][i]=(a[j][i]-sum((L[j][k]*L[i][k]*d[k] for k in range(i)),F(0)))/cur
    return True


class P92EMulticlusterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate=certify()

    def test_exact_scope_and_counts(self):
        c=self.certificate
        self.assertEqual(c['schema'],SCHEMA)
        self.assertEqual(c['case_count'],8)
        self.assertEqual(c['finite_exact_pair_count'],162)
        self.assertEqual(c['cluster_count_unconditional_residual_gate'],11)
        self.assertTrue(c['full_Preissmann_constant_for_up_to_11_uniform_clusters_under_gap_condition'])
        for field in ('all_arbitrary_spacings_global_Delta_bound_proved',
                      'arbitrary_number_of_clusters_global_Delta_bound_proved',
                      'full_Preissmann_Lemma6_proved','Riemann_Hypothesis_proved'):
            self.assertFalse(c[field])

    def test_pair_integral_upper_for_all_162(self):
        for case in self.certificate['cases']:
            for row in case['pairs']:
                finite=F(row['finite_HS_square'])
                upper=F(row['integral_HS_ceiling_square'])
                squpper=F(row['rational_sqrt_upper'])
                self.assertLessEqual(finite,upper)
                self.assertLessEqual(upper,F(1,8))
                self.assertGreaterEqual(squpper*squpper,upper)

    def test_uniform_pairs_zero(self):
        for s in (F(1,2),F(1),F(2)):
            self.assertEqual(two_cluster_coefficient(s,s),0)
            self.assertEqual(pair_hs_square_exact(s,s,4,10,s*s),0)
            self.assertEqual(pair_hs_integral_ceiling_square(s,s,s*s),0)

    def test_rational_sqrt_upper(self):
        for q in (F(0),F(1,8),F(13,17),F(1,10**80),F(101)):
            u=sqrt_rational_upper(q)
            self.assertGreaterEqual(u*u,q)
            if q>0:
                self.assertLessEqual((u-F(1,10**25))**2,q)

    def test_simple_real_point_cluster_model(self):
        ss=(F(1),F(1,2),F(3,2))
        sizes=(3,2,4)
        mult=(F(1),F(3,2))
        a,indices=exact_delta_matrix(ss,sizes,mult)
        self.assertTrue(all(a[i][j]==0 for i in range(len(a))
                            for j in range(len(a)) if indices[i]==indices[j]))
        self.assertTrue(all(a[i][j]>=0 for i in range(len(a)) for j in range(len(a))))
        for r in range(3):
            for s in range(r+1,3):
                hs_r=ss[r]**2; hs_s=ss[s]**2
                _,_,_,gr=build_uniform_clusters(ss,sizes,mult)
                gap=gr[s][0]-gr[r][-1]
                val=sum((a[i][j]**2 for i in range(len(a))
                         for j in range(len(a)) if indices[i]==r and indices[j]==s),F(0))
                self.assertEqual(val,pair_hs_square_exact(ss[r],ss[s],sizes[r],sizes[s],gap))

    def test_rational_ldlt_spectral_regression(self):
        for ss,sizes,mult in (
            ((F(1),F(1,2)),(2,3),(F(1),)),
            ((F(1),F(1,2),F(3,2)),(2,2,2),(F(1),F(1))),
            ((F(1),F(1,3),F(2),F(1,2)),(2,2,2,2),(F(2),F(3),F(1))),
        ):
            a,_=exact_delta_matrix(ss,sizes,mult)
            row=F(block_certificate(ss,sizes,mult)['row_sum_upper_exact'])
            bound=row+F(1,10**10)
            for sign in (-1,1):
                matrix=[[bound*(i==j)+sign*a[i][j] for j in range(len(a))]
                        for i in range(len(a))]
                self.assertTrue(exact_spd_ldlt(matrix))

    def test_up_to_11_clusters_and_conditional_14(self):
        self.assertGreater(F(25,8)**4,F(375,4))
        for row in self.certificate['cases']:
            if row['cluster_count']<=11:
                self.assertTrue(row['global_restricted_theorem_applies'])
                self.assertTrue(row['row_enclosure_below_target_sufficient_by_25_over_8'])
        fourteen=[x for x in self.certificate['cases'] if x['name']=='fourteen_far'][0]
        self.assertFalse(fourteen['global_restricted_theorem_applies'])
        self.assertTrue(fourteen['row_enclosure_below_target_sufficient_by_25_over_8'])

    def test_input_guards(self):
        with self.assertRaises(ValueError):
            build_uniform_clusters((F(1),F(2)),(1,2),(F(1),))
        with self.assertRaises(ValueError):
            build_uniform_clusters((F(1),F(2)),(2,2),(F(1,2),))
        with self.assertRaises(ValueError):
            pair_hs_integral_ceiling_square(F(1),F(3),F(1))
        with self.assertRaises(ValueError):
            sqrt_rational_upper(F(-1))

if __name__=='__main__':
    unittest.main()
