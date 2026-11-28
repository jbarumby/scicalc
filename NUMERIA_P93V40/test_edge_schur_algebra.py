import unittest

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from .edge_schur_algebra import TailLoadingError, finite_rank_schur


def _enclose(matrix, radius, *, skew_centers=False):
    entries = []
    for i in range(matrix.nrows()):
        for j in range(matrix.ncols()):
            shift = (1 if i < j else -1) * radius / 2 if skew_centers else 0
            entries.append(arb(matrix[i, j] + shift, radius))
    return arb_mat(matrix.nrows(), matrix.ncols(), entries)


def _eliminate_corrected_matrix(a, c, d, u, v, coefficient):
    """Build the corrected blocks and eliminate the tail in exact rationals."""
    shifted_head = a - coefficient * (u * u.transpose())
    shifted_cross = c - coefficient * (u * v.transpose())
    shifted_tail = d - coefficient * (v * v.transpose())
    elimination = shifted_tail.solve(shifted_cross.transpose())
    return shifted_head - shifted_cross * elimination, elimination


class FiniteRankSchurTests(unittest.TestCase):
    def test_rational_reference_inside_independently_rounded_symmetric_balls(self):
        with ctx.workprec(192):
            a = fmpq_mat([[4, fmpq(2, 7), -1], [fmpq(2, 7), 5, 1], [-1, 1, 7]])
            c = fmpq_mat([[fmpq(4, 11), -2], [1, fmpq(2, 13)], [3, 2]])
            d = fmpq_mat([[5, fmpq(-1, 3)], [fmpq(-1, 3), 4]])
            u = fmpq_mat([[2], [fmpq(-1, 5)], [3]])
            v = fmpq_mat([[1], [fmpq(2, 9)]])
            gamma = fmpq(1, 6)
            radius = fmpq(1, 2**24)
            enclosure = finite_rank_schur(
                _enclose(a, radius, skew_centers=True),
                _enclose(c, radius),
                _enclose(d, radius, skew_centers=True),
                _enclose(u, radius),
                _enclose(v, radius),
                arb(gamma, radius),
            )
            result = enclosure["schur"]
            for change in (fmpq(0), radius / 4, -radius / 4):
                perturbed_tail = d + fmpq_mat([[0, change], [change, 0]])
                perturbed_cross = c + fmpq_mat([[0, change], [0, 0], [0, 0]])
                expected, elimination = _eliminate_corrected_matrix(
                    a, perturbed_cross, perturbed_tail, u, v, gamma + change,
                )
                for i in range(expected.nrows()):
                    for j in range(expected.ncols()):
                        self.assertTrue(result[i, j].contains(expected[i, j]))
                        self.assertEqual(result[i, j].mid(), result[j, i].mid())
                        self.assertEqual(result[i, j].rad(), result[j, i].rad())
                for i in range(elimination.nrows()):
                    for j in range(elimination.ncols()):
                        self.assertTrue(enclosure["corrected_elimination"][i, j]
                                        .contains(elimination[i, j]))
            self.assertGreater(enclosure["denominator"], 0)
            self.assertGreater(enclosure["tail_loading"], 0)

    def test_endpoint_in_coupling_range_cancels_the_schur_update(self):
        a, c, d = arb_mat([[7]]), arb_mat([[6, 5]]), arb_mat([[2, 0], [0, 4]])
        u, v = arb_mat([[8]]), arb_mat([[1], [4]])
        enclosure = finite_rank_schur(a, c, d, u, v, arb(1) / 8)
        self.assertTrue(enclosure["effective_endpoint"][0, 0].is_zero())
        self.assertEqual(enclosure["base_schur"][0, 0], enclosure["schur"][0, 0])
        self.assertEqual(enclosure["tail_loading"], arb(9) / 16)

    def test_singular_negative_or_indeterminate_corrected_tail_is_rejected(self):
        a = d = v = arb_mat([[1]])
        c = u = arb_mat([[0]])
        for gamma in (arb(1), arb(2), arb(1, "1e-8")):
            with self.subTest(gamma=str(gamma)), self.assertRaises(TailLoadingError) as failure:
                finite_rank_schur(a, c, d, u, v, gamma)
            self.assertFalse(failure.exception.denominator > 0)
        with self.assertRaisesRegex(ArithmeticError, "solve"):
            finite_rank_schur(a, c, arb_mat([[0]]), u, v, arb(0))

    def test_incompatible_inputs_are_rejected(self):
        a = d = u = v = arb_mat([[1]])
        c = arb_mat([[0]])
        with self.assertRaisesRegex(ValueError, "coefficient"):
            finite_rank_schur(a, c, d, u, v, arb(-1))
        with self.assertRaisesRegex(ValueError, "dimensions"):
            finite_rank_schur(a, arb_mat([[0, 0]]), d, u, v, arb(0))
        with self.assertRaises(ValueError):
            finite_rank_schur(a, arb_mat([[0, 0]]), arb_mat([[1, 2], [0, 1]]),
                              u, arb_mat([[1], [1]]), arb(0))


if __name__ == "__main__":
    unittest.main()
