"""Independent strict-scope tests for P92R exact arithmetic (not a formal proof)."""
import unittest
from fractions import Fraction as F
from .verify import SCHEMA,verify,lemma5,harmonic,ramp_U0,ramp_bellman,l6,finite_compensator,nearest,pi_lower


class P92RTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=verify()

    def test_census_and_scientific_guards(self):
        p=self.cert
        self.assertEqual(p["schema"],SCHEMA)
        self.assertEqual(p["counts"],{"L5":222,"ramp":400,"mixed":12,"L6":60})
        self.assertTrue(p["analytic_independent_L5_proof_in_README"])
        self.assertTrue(p["analytic_general_L6_reconstruction_proposed_for_audit"])
        self.assertFalse(p["independent_peer_review_completed"])
        self.assertFalse(p["Riemann_Hypothesis_proved"])

    def test_uniform_harmonic_identity(self):
        for n in range(1,35):
            p=lemma5([F(1)]*n)
            self.assertEqual(F(p["sum"]),harmonic(n))
            self.assertEqual(F(p["margin"]),0)

    def test_integer_endpoint_nearly_sharp_remainder(self):
        for eps in (F(1,100),F(1,1000),F(1,10000)):
            p=lemma5([F(1),F(1)-eps])
            self.assertEqual(p["harmonic_index"],1)
            self.assertGreater(F(p["sum"]),1)
            self.assertLess(F(p["margin"]),3*eps)

    def test_ramp_formula_against_direct_enumeration(self):
        for h in (F(1,3),F(1),F(5,4)):
            for L in (F(0),F(1,3),F(2),F(7,2)):
                n=int(L//h)+int(L%h>0)
                exact=sum((h*max(F(0),L-j*h) for j in range(n+1)),F(0))
                self.assertEqual(exact,ramp_U0(L,h))

    def test_bellman_both_branches(self):
        for L in (F(0),F(1),F(8,3),F(12)):
            for h in (F(1,3),F(2),F(4)):
                for b in (F(1,5),F(1),F(3),F(10)):
                    self.assertGreaterEqual(ramp_bellman(L,h,b),0)

    def test_nearest_neighbor_endpoint_convention(self):
        z,d=nearest([F(1),F(2),F(1,3),F(5,2)])
        self.assertEqual(z,[F(1),F(3),F(10,3),F(35,6)])
        self.assertEqual(d,[F(1),F(1,3),F(1,3),F(5,2)])

    def test_mixed_cutoff_remainder_decay(self):
        g=[F(1),F(3),F(1,2),F(2)]
        a=finite_compensator(g,10,F(1),F(1,2))
        b=finite_compensator(g,100,F(1),F(1,2))
        self.assertLess(F(b["remainder"]),F(a["remainder"]))
        self.assertGreaterEqual(F(b["margin"]),0)

    def test_L6_affine_invariance(self):
        p=[F(0),F(1),F(9,5),F(2),F(11)]
        a=l6(p,0,len(p)-1)
        scale=F(7,13)
        b=l6([scale*z+F(37,10) for z in p],0,len(p)-1)
        self.assertEqual(F(b["lhs"]),F(a["lhs"])/scale**3)
        self.assertEqual(F(b["rhs_pi_lower"]),F(a["rhs_pi_lower"])/scale**3)

    def test_L6_interior_pair(self):
        p=[F(0),F(2),F(5,2),F(21,5),F(7),F(11)]
        a=l6(p,1,len(p)-2)
        b=l6([-z for z in reversed(p)],1,len(p)-2)
        self.assertEqual(F(a["lhs"]),F(b["lhs"]))
        self.assertGreaterEqual(F(a["margin"]),0)

    def test_all_exportable_margins_and_pi_lower(self):
        self.assertGreater(pi_lower(),3)
        self.assertLess(pi_lower(),F(22,7))
        for p in self.cert["l5"]:
            self.assertGreaterEqual(F(p["margin"]),0)
        for p in self.cert["bellman"]:
            self.assertGreaterEqual(F(p),0)
        for p in self.cert["mixed"]:
            self.assertGreaterEqual(F(p["margin"]),0)
        for p in self.cert["l6"]:
            self.assertGreaterEqual(F(p["margin"]),0)


if __name__ == "__main__":
    unittest.main()
