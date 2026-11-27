"""Independent regression tests of the block interval LDL."""

from fractions import Fraction
import unittest

from flint import arb, arb_mat, ctx

from .interval_ldl import certify_inertia, certify_prefixes


def _fraction(dyadic):
    m, e = map(int, dyadic)
    return m * Fraction(2)**e


def _reference_pivots(entries):
    """Exact determinant ratios by scalar elimination, independent of blocks."""
    matrix = [list(map(Fraction, row)) for row in entries]
    pivots = []
    for k, row in enumerate(matrix):
        pivot = row[k]
        pivots.append(pivot)
        for i in range(k + 1, len(matrix)):
            for j in range(k + 1, len(matrix)):
                matrix[i][j] -= matrix[i][k] * row[j] / pivot
    return pivots


class PrefixLDLTests(unittest.TestCase):
    def test_dense_exact_crosscheck_and_checkpoint_reuse(self):
        # Irregular exact rational L and signed D, with a nontrivial trailing update.
        lower = [[1, 0, 0, 0, 0], [3, 1, 0, 0, 0], [2, -1, 1, 0, 0],
                 [1, 2, -3, 1, 0], [7, 2, 1, -4, 1]]
        weights = [2, 3, 5, 7, 11]
        entries = [[sum(lower[i][k] * weights[k] * lower[j][k] for k in range(5))
                    for j in range(5)] for i in range(5)]
        exact = _reference_pivots(entries)
        with ctx.workprec(256):
            for block in (1, 2, 3, 8):
                proof = certify_prefixes(arb_mat(entries), (1, 3, 5), block_size=block)
                self.assertTrue(proof["certified"])
                self.assertEqual(proof["pivots_proved"], 5)
                for target, (lo, hi) in zip(exact, proof["pivot_intervals_dyadic"], strict=True):
                    self.assertLessEqual(_fraction(lo), target)
                    self.assertGreaterEqual(_fraction(hi), target)

    def test_negative_shifted_direction_is_not_negative_unshifted(self):
        result = certify_prefixes(arb_mat([[2, 0], [0, 1]]), (1, 2), margin="3/2")
        self.assertEqual(result["first_obstruction"]["pivot_index"], 2)
        self.assertEqual(result["prefixes"][0]["status"], "certified")
        self.assertEqual(result["prefixes"][1]["status"], "negative_shifted_direction")
        self.assertFalse(result["prefixes"][1]["negative_direction_of_unshifted_matrix"])

    def test_indeterminate_and_strict_negative_are_distinct(self):
        uncertain = certify_prefixes(arb_mat([[1, 0], [0, arb(0, 1)]]), (1, 2))
        self.assertEqual(uncertain["status"], "indeterminate_pivot")
        self.assertEqual(uncertain["prefixes"][1]["status"], "inconclusive")
        negative = certify_prefixes(arb_mat([[1, 0], [0, -1]]), (1, 2))
        self.assertEqual(negative["status"], "unexpected_pivot_sign")
        self.assertTrue(negative["prefixes"][1]["negative_direction_of_unshifted_matrix"])

    def test_inertia_has_both_signs(self):
        matrix = arb_mat([[-3, 1, 2], [1, -4, 1], [2, 1, 2]])
        proof = certify_inertia(matrix, (-1, -1, 1), block_size=2)
        self.assertTrue(proof["certified"])
        self.assertEqual(proof["inertia"], {"positive": 1, "negative": 2, "zero": 0})

    def test_reject_bad_input(self):
        symmetric = arb_mat([[1, 0], [0, 1]])
        with self.assertRaisesRegex(ValueError, "asymmetric"):
            certify_prefixes(arb_mat([[1, 0], [1, 1]]), (2,))
        with self.assertRaises(ValueError):
            certify_prefixes(symmetric, (2, 1))
        with self.assertRaises(ValueError):
            certify_prefixes(symmetric, (1,), 1e-60)
        with self.assertRaises(ValueError):
            certify_inertia(symmetric, (-1, True))


if __name__ == "__main__":
    unittest.main()
