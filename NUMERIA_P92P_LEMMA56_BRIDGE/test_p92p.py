"""P92P: regression tests. General theorems are in README, not in these samples."""
import unittest
from fractions import Fraction as F
from .certificate import audit, certify, pi_bounds, gaps, energy, potential, ldlt

class P92PTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=certify()
        cls.c=F(cls.result['pi_lower'])**2/3

    def test_scope_is_honest(self):
        d=self.result
        self.assertTrue(d['universal_algebraic_bridge_proved'])
        self.assertTrue(d['rational_schur_equivalence_proved'])
        for field in ('all_configurations_lemma5_proved',
                      'all_configurations_lemma6_independently_proved',
                      'Riemann_Hypothesis_proved'):
            self.assertFalse(d[field])

    def test_pi_lower_machin(self):
        a,b=pi_bounds()
        self.assertLess(F(3),a)
        self.assertLess(a,b)
        self.assertLess(b,F(22,7))

    def test_63_fraction_cases_and_27_ldl(self):
        d=self.result
        self.assertEqual(d['exact_fraction_case_count'],63)
        self.assertEqual(d['independent_dense_ldlt_count'],27)
        self.assertEqual(len(d['cases']),63)
        self.assertEqual(sum(c['ldlt_min_pivot'] is not None for c in d['cases']),27)
        for c in d['cases']:
            self.assertGreater(F(c['slack']),0)
            self.assertEqual(F(c['rhs_budget'])-F(c['h_squared_kernel']),F(c['slack']))
            if c['ldlt_min_pivot'] is not None:
                self.assertGreater(F(c['ldlt_min_pivot']),0)

    def test_compensation_identity_by_direct_sums(self):
        xx=[F(t) for t in [0,2,7,10,100,101]]
        d=gaps(xx)
        s,t=0,4
        h=xx[t]-xx[s]
        direct=sum((d[i]*h*h/((xx[i]-xx[s])**2*(xx[i]-xx[t])**2)
                    for i in range(len(xx)) if i not in (s,t)),F(0))
        other=energy(xx,d,s)+energy(xx,d,t)-2*(potential(xx,d,s)-potential(xx,d,t))/h
        self.assertEqual(direct+3*(d[s]+d[t])/h**2,other)

    def test_index_separation_five_sufficient_sign(self):
        z=audit([0,7,8,16,17,32,33],0,5,self.c,True)
        self.assertTrue(z['monotone_potential_suffices'])
        self.assertEqual(F(z['potential_difference']),F(31333,261800))
        self.assertEqual(z['index_separation'],5)

    def test_index_separation_six_sufficient_sign(self):
        z=audit([0,8000,8300,32300,32500,32510,33510,33515],0,6,self.c,True)
        self.assertGreater(F(z['potential_difference']),F(47,100))
        self.assertEqual(z['index_separation'],6)

    def test_monotone_potentials_not_necessary(self):
        z=audit(list(range(13)),2,9,self.c,True)
        self.assertFalse(z['monotone_potential_suffices'])
        self.assertGreater(F(z['slack']),0)

    def test_ldl_detects_false_schur_bound(self):
        m=[[F(1),F(0),F(1)],[F(0),F(1),F(1)],[F(1),F(1),F(1)]]
        self.assertLess(ldlt(m)[-1],0)

if __name__=='__main__':
    unittest.main()
