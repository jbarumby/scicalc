"""P92Q: independent regressions for triangular packing variational lemma."""
import unittest
from fractions import Fraction as F
from .certificate import (
    SCHEMA,certify,exact_case,interior_nearest_gaps,tent,mesh_kernel,
    D_mesh,harmonic,pi_lower,full_L6_left_right
)

class P92QTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_census_and_claim_scope(self):
        d=self.data
        self.assertEqual(d["schema"],SCHEMA)
        self.assertEqual(d["exact_rational_regressions"],85)
        self.assertEqual(d["full_L6_checks_with_rational_pi_lower"],5)
        self.assertTrue(d["analytic_universal_for_stated_lattice_class"])
        self.assertFalse(d["general_Lemma5_and_6_for_arbitrary_real_positions"])
        self.assertFalse(d["general_Riemann_Hypothesis_proved"])

    def test_tents_mass_centre_occupancy(self):
        for M,J in ((9,(2,5)),(20,(3,11,18)),(100,(1,2,5,19,20,46,73,74,99))):
            g=interior_nearest_gaps(M,J)
            for j,d in g.items():
                self.assertEqual(sum((tent(j,d,k) for k in range(1,M)),F(0)),d)
                self.assertEqual(sum((k*tent(j,d,k) for k in range(1,M)),F(0)),j*d)
            for k in range(1,M):
                self.assertLessEqual(sum((tent(j,d,k) for j,d in g.items()),F(0)),1)

    def test_sparse_convex_kernel_majorization(self):
        for M,J in ((7,(1,3,6)),(15,(3,7,12)),(100,(1,2,5,19,20,46,73,74,99))):
            g=interior_nearest_gaps(M,J)
            lhs=sum((d*mesh_kernel(M,j) for j,d in g.items()),F(0))
            rhs=sum((mesh_kernel(M,k) for k in range(1,M)),F(0))
            self.assertLessEqual(lhs,rhs)

    def test_general_convex_quadratics(self):
        for M,J in ((6,(1,4)),(10,(3,7)),(20,(3,10,15))):
            g=interior_nearest_gaps(M,J)
            for c in (F(0),F(M,2),F(M,3),F(M)):
                f=lambda k:(F(k)-c)**2
                self.assertLessEqual(sum((d*f(j) for j,d in g.items()),F(0)),
                                     sum((f(k) for k in range(1,M)),F(0)))

    def test_harmonic_identity_in_every_degree(self):
        for M in range(1,55):
            partial=sum((mesh_kernel(M,k) for k in range(1,M)),F(0))
            self.assertEqual(M**3*partial,2*D_mesh(M))
            self.assertEqual(M**3*partial,2*M*harmonic(M-1,2)+4*harmonic(M-1))

    def test_thinning_slack_is_positive_for_irregular(self):
        t=exact_case(100,(1,2,5,19,20,46,73,74,99))
        self.assertGreater(F(t["thinning_margin_exact"]),0)
        self.assertEqual(F(t["normalized_interior_sum_exact"])+
                         F(t["thinning_margin_exact"]),F(t["normalized_complete_sum_exact"]))

    def test_full_sharp_L6_with_outside_points(self):
        for M,J in ((8,(2,3,6)),(17,(1,6,10,16)),(40,(1,2,20,39))):
            x=[F(-1,2),F(0)]+[F(j) for j in J]+[F(M),F(M+1,3)]
            lhs,rhs,d=full_L6_left_right(x,1,len(x)-2)
            self.assertLess(lhs,rhs)
            self.assertLessEqual(d[1],1)
            self.assertLessEqual(d[-2],1)

    def test_affine_invariance(self):
        x=[F(-1),F(0),F(1),F(2),F(5),F(19),F(20),F(46),F(73),F(74),F(99),F(100),F(101)]
        s=1;t=len(x)-2
        lo,ro,_=full_L6_left_right(x,s,t)
        scale=F(3,7)
        shift=F(41,11)
        l1,r1,_=full_L6_left_right([scale*z+shift for z in x],s,t)
        self.assertEqual(l1,lo/scale**3)
        self.assertEqual(r1,ro/scale**3)

    def test_critical_input_validation(self):
        with self.assertRaises(ValueError): interior_nearest_gaps(10,(5,5))
        with self.assertRaises(ValueError): interior_nearest_gaps(10,(0,3))
        with self.assertRaises(ValueError): interior_nearest_gaps(0,())
        with self.assertRaises(ValueError): mesh_kernel(9,9)
        with self.assertRaises(ValueError): tent(3,0,4)

    def test_machin_lower_and_certified_scope(self):
        self.assertGreater(pi_lower(),3)
        self.assertLess(pi_lower(),F(22,7))
        for c in self.data["cases"]:
            self.assertGreaterEqual(F(c["thinning_margin_exact"]),0)
            self.assertEqual(F(c["twice_D_of_1_over_M_exact"]),
                             F(c["normalized_complete_sum_exact"]))
        for c in self.data["full_checks"]:
            self.assertGreater(F(c["full_L6_margin_with_pi_lower_exact"]),0)


if __name__=="__main__":
    unittest.main()
