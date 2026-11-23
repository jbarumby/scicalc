"""P92-M exact Fraction regressions; the general proof is in README.md."""
from fractions import Fraction as F
import unittest
from .certificate import (
    SCHEMA, T, audit_case, make_configuration, delta_matrix, exact_ldl_positive,
    obstruction, certify
)

class P92MTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_scope_and_census(self):
        x=self.data
        self.assertEqual(x['schema'],SCHEMA)
        self.assertEqual(x['rational_configurations'],10)
        self.assertEqual(x['full_exact_ldl'],7)
        self.assertEqual(len(x['sharp_exponent_obstructions']),6)
        self.assertGreater(x['rational_pair_checks'],100)
        for guard in ('unrestricted_point_configurations_proved','general_Preissmann_Lemma6_proved','Riemann_Hypothesis_proved'):
            self.assertFalse(x[guard],guard)

    def test_rational_Machin_comparison(self):
        p=F(self.data['exact_machin_pi_lower'])
        self.assertGreater(p,F(25,8))
        self.assertGreater(F(25,8)**4,F(375,4))
        self.assertLess(F(25,2),F(2,15)*p**4)

    def test_geometry_and_distances(self):
        configs=(((F(2),F(1,3)),(2,7)),((F(1),F(1,2),F(5)),(6,3,2)),((F(1),),(9,)))
        for roots,sizes in configs:
            xs, ids, group, rr, steps = make_configuration(roots,sizes,3)
            self.assertEqual(len(xs),sum(sizes))
            for i,x in enumerate(xs):
                d=min(abs(x-y) for j,y in enumerate(xs) if i!=j)
                self.assertEqual(d,steps[ids[i]])

    def test_zero_interactions_equal_steps(self):
        xs,ids,groups,roots,steps=make_configuration([1,1,1],[2,3,4])
        for row in delta_matrix(xs,ids,roots):
            self.assertTrue(all(v==0 for v in row))

    def test_sample_two_and_three_blocks(self):
        a=audit_case('mini',[F(1),F(1,6)],[7,9],3,True)
        self.assertEqual(a['K'],2)
        self.assertEqual(a['pair_interactions'],1)
        self.assertGreater(F(a['comparison_ldl_min_pivot']),0)
        self.assertGreater(F(a['full_ldl_min_pivot']),0)
        b=audit_case('triple',[F(1),F(1,3),F(3)],[3,5,4],3,True)
        self.assertEqual(b['pair_interactions'],3)
        self.assertGreater(F(b['full_ldl_min_pivot']),0)

    def test_all_case_LDL_and_pair_checks(self):
        for case in self.data['cases']:
            self.assertGreater(F(case['comparison_ldl_min_pivot']),0)
            if case['full_ldl_min_pivot'] is not None:
                self.assertGreater(F(case['full_ldl_min_pivot']),0)
            self.assertEqual(case['pair_interactions'],case['K']*(case['K']-1)//2)

    def test_sharp_exponent_obstruction(self):
        for case in self.data['sharp_exponent_obstructions']:
            n=case['n']
            lower=F(case['analytic_lower_on_operator_norm'])
            self.assertEqual(lower,(1-F(1,n))**2/(18*n))
            self.assertGreaterEqual(lower,F(1,72*n))
            self.assertEqual(case['rho'],n)
            self.assertEqual(case['m_fine'],n**3)
            if case['bilinear_exact'] is not None:
                self.assertGreater(F(case['bilinear_exact']),lower)

    def test_independent_obstruction_bilinear(self):
        for n in (2,3,4):
            case=obstruction(n,exact_sum=True)
            self.assertGreater(F(case['bilinear_exact']),F(case['analytic_lower_on_operator_norm']))

    def test_ldl_known_positive_matrix(self):
        M=[[F(2),F(1,2),F(0)],[F(1,2),F(1),F(1,4)],[F(0),F(1,4),F(3)]]
        self.assertGreater(exact_ldl_positive(M),0)
        with self.assertRaises(ArithmeticError):
            exact_ldl_positive([[0,F(1)],[F(1),0]])

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            make_configuration([F(0)],[2])
        with self.assertRaises(ValueError):
            make_configuration([1],[1])
        with self.assertRaises(ValueError):
            make_configuration([1,2],[3])
        with self.assertRaises(ValueError):
            obstruction(1)

if __name__=='__main__':
    unittest.main()
