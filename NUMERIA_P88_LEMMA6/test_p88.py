"""P88 exact Fraction tests; source-uniform claims are proven in README."""
import unittest
from fractions import Fraction as F
from .compensation import (
    SCHEMA, reciprocal_kernel_identity, harmonics, uniform_cross_parts,
    lattice_finite, mixed_balance, audit, machin_pi_lower,
)

class P88Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc=audit()

    def test_lemmas_not_promoted_falsely(self):
        d=self.doc
        self.assertEqual(d["schema"], SCHEMA)
        self.assertTrue(d["uniform_sharp_coefficient_negative_three"])
        self.assertTrue(d["universal_relaxed_mixed_inequality_proved"])
        self.assertFalse(d["arbitrary_spacing_sharp_compensation_proved"])
        self.assertFalse(d["independent_full_Preissmann_Lemma6_proved"])
        self.assertFalse(d["RZ_S_epsilon_reproved"])
        self.assertFalse(d["RH_proved"])

    def test_partial_fraction_identity_exact(self):
        for u,v in ((F(1,7),-F(13,11)),(-F(1,2),F(3,7)),(F(100),F(1,200))):
            c=reciprocal_kernel_identity(u,v)
            self.assertEqual(F(c["lhs"]),F(c["rhs"]))

    def test_cross_harmonic_identity_for_any_integer_m_samples(self):
        for m in (1,2,3,4,5,8,13,20,40):
            r=uniform_cross_parts(m)
            self.assertEqual(F(r["cross_sum_total"]),F(2,m*m))
            self.assertEqual(F(r["left_tail"]),harmonics(m)/m)
            self.assertEqual(F(r["between_endpoints"]),-2*harmonics(m-1)/m)

    def test_lattice_tail_is_rigorously_bounded(self):
        self.assertEqual(len(self.doc["finite_lattice_tail_examples"]),7)
        for row in self.doc["finite_lattice_tail_examples"]:
            self.assertTrue(row["exact_tail_bound_verified"])
            lo,hi=map(F,(row["limiting_minus_finite_lower"],row["limiting_minus_finite_upper"]))
            self.assertTrue(0<lo<hi<=F(row["positive_tail_upper_exact"]))

    def test_strong_Finite_checks_all_perturbed_configurations(self):
        self.assertGreater(self.doc["finite_mixed_pair_count"],700)
        self.assertEqual(sum(v["pairs"] for v in self.doc["finite_mixed_configurations"]),
                         self.doc["finite_mixed_pair_count"])
        for c in self.doc["finite_mixed_configurations"]:
            self.assertGreater(F(c["minimum_positive_compensation_with_pi_lower"]),0)
            for row in c["records"]:
                h=F(row["endpoint_distance_exact"])
                E=F(row["compensation_with_pi_lower_exact"])
                R=F(row["sharp_RHS_with_rational_pi_lower"])
                S=F(row["sharp_original_LHS_exact"])
                self.assertEqual(E,h*h*(R-S))
                self.assertGreater(E,0)
                self.assertGreaterEqual(F(row["relaxed_RHS_with_pi_lower"]),S)

    def test_counterexample_cross_can_be_negative(self):
        row=mixed_balance(["0","1","2","3"],0,3)
        self.assertEqual(F(row["mixed_cross_exact"]),-1)
        self.assertTrue(F(row["compensation_with_pi_lower_exact"])>0)

    def test_Machin_lower_is_not_float(self):
        pi=machin_pi_lower()
        self.assertGreater(pi,F(157,50))
        self.assertLess(pi,F(22,7))

    def test_reject_degenerate_and_impossible_input(self):
        with self.assertRaises(ValueError):
            reciprocal_kernel_identity(1,1)
        with self.assertRaises(ValueError):
            reciprocal_kernel_identity(0,1)
        with self.assertRaises(ValueError):
            lattice_finite(0,2)
        with self.assertRaises(ValueError):
            lattice_finite(2,0)
        with self.assertRaises(ValueError):
            lattice_finite(2,4,0)
        with self.assertRaises(ValueError):
            uniform_cross_parts(0)
        with self.assertRaises(ValueError):
            mixed_balance(["0","1","2"],0,0)
        with self.assertRaises(ValueError):
            mixed_balance(["0","1","2"],0,3)

if __name__=="__main__":
    unittest.main()
