"""P256: exact rational and fail-closed regression tests."""
import unittest
from fractions import Fraction as Q
from .weighted_bridge import (
    prime_power_base, prime_terms, weighted_cosine_coefficients,
    abel_coefficients, imaginary_mellin_coefficients, grouped_log_coefficients,
    slope_coefficients, activation_7, git_blob_sha1, audit, validate_source_tree)

class WeightedTests(unittest.TestCase):
    def test_prime_powers(self):
        self.assertEqual([prime_power_base(n) for n in range(2,13)],
                         [2,3,2,5,None,7,2,3,None,11,None])

    def test_four_term_window(self):
        self.assertEqual(prime_terms(Q(6)),[(2,2),(3,3),(4,2),(5,5)])
        self.assertEqual(grouped_log_coefficients(Q(6)),
                         {2:Q(5,2),3:Q(1),5:Q(1,5)})

    def test_strict_boundary(self):
        self.assertNotIn(7,weighted_cosine_coefficients(Q(7)))
        self.assertIn(7,weighted_cosine_coefficients(Q(7)+Q(1,10000)))
        self.assertNotIn(8,weighted_cosine_coefficients(Q(8)))
        self.assertIn(8,weighted_cosine_coefficients(Q(8)+Q(1,2)))

    def test_abel_coefficient_identity(self):
        for x in [Q(2),Q(7),Q(71,10),Q(10),Q(999,10),Q(100)]:
            with self.subTest(x=x):
                self.assertEqual(weighted_cosine_coefficients(x),abel_coefficients(x))

    def test_imaginary_mellin_identity(self):
        for x in [Q(3),Q(6),Q(7),Q(701,100),Q(22),Q(77)]:
            with self.subTest(x=x):
                self.assertEqual(weighted_cosine_coefficients(x),imaginary_mellin_coefficients(x))

    def test_7_new_term_small(self):
        for eps in [Q(1,10),Q(1,10000),Q(1,1000000)]:
            with self.subTest(eps=eps):
                self.assertEqual(activation_7(eps)['extra_weighted_log7_cos_coefficient'],str(eps/7))

    def test_increase_is_linear_at_7(self):
        eps=Q(1,100)
        weighted=weighted_cosine_coefficients(Q(7)+eps)
        self.assertEqual(weighted[7],(7,eps/7))
        self.assertEqual(set(weighted),{2,3,4,5,7})

    def test_slope_jump_7(self):
        left=slope_coefficients(Q(7))
        right=slope_coefficients(Q(7)+Q(1,100))
        self.assertEqual(right[7],Q(1,7))
        for p in left:
            self.assertEqual(left[p],right[p])

    def test_slope_jump_4(self):
        left=slope_coefficients(Q(4))
        right=slope_coefficients(Q(4)+Q(1,100))
        self.assertEqual(right[2]-left[2],Q(1,4))

    def test_nonprime_6_no_jump(self):
        self.assertEqual(slope_coefficients(Q(6)),slope_coefficients(Q(6)+Q(1,100)))

    def test_zero_and_negative_unsupported(self):
        for x in [Q(0),Q(-1)]:
            with self.assertRaises(ValueError):
                weighted_cosine_coefficients(x)

    def test_invalid_prime_input(self):
        for n in [1,0,-1,7.0]:
            with self.assertRaises(ValueError):
                prime_power_base(n)

    def test_epsilon_out_of_range(self):
        for eps in [0,-1,1,Q(3,2)]:
            with self.assertRaises(ValueError):
                activation_7(Q(eps))

    def test_summary_does_not_claim_global_rh(self):
        r=audit()
        self.assertFalse(r['RH_proved'])
        self.assertFalse(r['new_Weil_coercivity_proved'])
        self.assertEqual(len(r['tests_in_exact_fraction_arithmetic']),9)

    def test_blob_hash_changes_under_mutation(self):
        self.assertNotEqual(git_blob_sha1(b'abc'),git_blob_sha1(b'abd'))

    def test_source_absence_fails_closed(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        with TemporaryDirectory() as t:
            with self.assertRaises(FileNotFoundError):
                validate_source_tree(Path(t))

if __name__=='__main__':
    unittest.main()
