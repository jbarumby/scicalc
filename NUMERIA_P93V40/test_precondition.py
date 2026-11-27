"""Exercise exact congruence and distinguish shifted from unshifted failures."""

from fractions import Fraction
import unittest

from flint import arb, arb_mat, ctx

from .precondition import certify_by_congruence, verify_triangular_congruence


class PreconditionTests(unittest.TestCase):
    def test_high_conditioning(self):
        with ctx.workprec(384):
            epsilon = arb(2)**-200
            a = arb_mat([[1, 1], [1, 1 + epsilon]])
            result = certify_by_congruence(a, [1,2], Fraction(1, 2**205))
            self.assertEqual(result["status"], "certified")
            self.assertTrue(all(row["certified"] for row in result["prefixes"]))

    def test_shifted_failure_has_exact_witness(self):
        with ctx.workprec(256):
            a = arb_mat([[1,0], [0,1]])
            result = certify_by_congruence(a, [1,2], Fraction(3,2))
            self.assertEqual(result["status"], "strictly_below_requested_margin")
            self.assertFalse(result["rayleigh_check"]["proven_negative"])
            self.assertTrue(result["rayleigh_check"]["witness_dyadic"])

    def test_reject_nontriangular_proposal(self):
        with self.assertRaises(ValueError):
            verify_triangular_congruence(arb_mat([[1,0],[0,1]]),
                                         arb_mat([[1,0],[1,1]]), [2], 0)
