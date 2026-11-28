"""Tests for the analytic infinite budget and the rigorous finite margin."""

from fractions import Fraction
import unittest

from flint import arb, arb_mat, ctx

from .edge_schur import _proof, complement_budget, edge_constants


def _ratio(encoded):
    mantissa, exponent = encoded
    return int(mantissa) * Fraction(2)**exponent


class EdgeSchurTests(unittest.TestCase):
    def test_budget_prices_rank_and_scalar_remainder_separately(self):
        with ctx.workprec(512):
            rank, remainder, slack = edge_constants()
            self.assertTrue(rank > 0)
            self.assertTrue(0 < remainder < arb(2)**-680)
            self.assertTrue(slack > arb(1) / 200)

    def test_absolute_margin_is_in_the_orthonormal_basis(self):
        with ctx.workprec(256):
            proof = _proof(arb_mat([[2, 1], [1, 2]]))
            self.assertEqual(proof["status"], "certified")
            # The eigenvalues of this test matrix are exactly 1 and 3.
            self.assertTrue(0 < _ratio(proof["absolute_margin_lower_dyadic"]) <= 1)

    def test_infinite_reserve_and_calculated_head_requirements(self):
        report = complement_budget(384)
        rows = {(row["modes"], row["parity"]): row for row in report["results"]}
        for parity in (0, 1):
            self.assertFalse(rows[400, parity]["tail_certified_by_bound"])
            large = rows[3309, parity]
            self.assertTrue(large["tail_certified_by_bound"])
            self.assertGreater(_ratio(large["corrected_tail_reserve_lower_dyadic"]),
                               Fraction("0.118"))
            later = rows[3400, parity]
            self.assertLess(_ratio(later["required_corrected_head_gap_upper_dyadic"]),
                            _ratio(large["required_corrected_head_gap_upper_dyadic"]))


if __name__ == "__main__":
    unittest.main()
