"""Independent finite tests for the V.39 harmonic congruence in Arb."""

from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from flint import arb, arb_mat, ctx

from .harmonic import THETAS, assemble_congruence, certify_hybrid, conditional_schur


def b(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def matrix(rows):
    return arb_mat([[b(x) for x in row] for row in rows])


class HarmonicTests(unittest.TestCase):
    def test_five_congruences_against_independent_exact_inverse(self):
        with ctx.workprec(192):
            a = matrix([[2, F(1, 5), F(1, 6)], [F(1, 5), 3, F(1, 7)], [F(1, 6), F(1, 7), 4]])
            penalty = matrix([[F(1, 10)]])
            hq = matrix([[F(1, 3), F(1, 7)], [F(1, 7), F(2, 3)]])
            c, d, lam, margin = F(19, 10), F(21, 5), F(91, 100), F(1, 10**60)
            s = lam * d / ((1 - lam) * c)
            r, t, v = 1 + s / 3, s / 7, 1 + 2 * s / 3
            j = matrix([[1, 0, 0], [0, r / (1 + s), t / (1 + s)],
                        [0, t / (1 + s), v / (1 + s)]])
            det = r * v - t**2
            inverse = matrix([[v / det, -t / det], [-t / det, r / det]])
            f = (matrix([[1, 0], [0, 1]]) - inverse) * b(c / lam)

            for theta in THETAS:
                with self.subTest(theta=theta):
                    expected_head = arb_mat(a)
                    expected_head[0, 0] -= b(F(1, 10))
                    for i in range(3):
                        expected_head[i, i] -= b(margin)
                    for i in range(2):
                        expected_head[i + 1, i + 1] -= b((1 - theta) * c)
                        for k in range(2):
                            expected_head[i + 1, k + 1] -= b(theta) * f[i, k]
                    expected = j * expected_head * j
                    actual = assemble_congruence(a, penalty, hq, c, d, lam, theta, margin)
                    for i in range(3):
                        for k in range(3):
                            self.assertTrue(actual[i, k].overlaps(expected[i, k]))

    def test_localized_prefix_shares_one_computation(self):
        with ctx.workprec(192):
            a, penalty = matrix([[2, F(1, 3)], [F(1, 3), 2]]), matrix([[1, 0], [0, 1]])
            report = certify_hybrid(a, penalty, None, "19/10", "21/5", "91/100", "1e-60")
            self.assertEqual(report["distinct_assembled_heads"], 1)
            self.assertEqual([r["theta_exact"] for r in report["results"]], [str(t) for t in THETAS])
            first = report["results"][0]["ldl"]
            self.assertTrue(first["certified"])
            self.assertTrue(all(r["ldl"] == first for r in report["results"]))

    def test_hybrid_changes_result_on_a_trailing_mode(self):
        with ctx.workprec(192):
            a, penalty, hq = matrix([[2, 0], [0, 1]]), matrix([[F(1, 10)]]), matrix([[F(1, 100)]])
            report = certify_hybrid(a, penalty, hq, "19/10", "21/5", "91/100", "1e-60")
            self.assertEqual(report["distinct_assembled_heads"], 5)
            self.assertEqual([r["ldl"]["certified"] for r in report["results"]],
                             [False, False, False, True, True])

    def test_reject_invalid_congruence_parameters(self):
        a, penalty, hq = matrix([[2, 0], [0, 1]]), matrix([[1]]), matrix([[1]])
        with self.assertRaises(ValueError):
            assemble_congruence(a, penalty, hq, 1, 2, "91/100", 0.5)
        with self.assertRaises(ValueError):
            assemble_congruence(a, penalty, None, 1, 2, "91/100", "1/2")

    def test_compare_conditional_schur_thresholds(self):
        source = Path(__file__).resolve().parents[1] / "NUMERIA_P93V39"
        parameters = json.loads((source / "hybrid39.json").read_text(encoding="utf-8"))
        with ctx.workprec(192):
            report = conditional_schur(parameters["c"], parameters["d"], parameters["lambda"])
        self.assertTrue(report["conditional_on_full_shifted_heads"])
        rows = report["results"]
        self.assertEqual([rows[i]["first_sufficient_modes_for_this_budget"] for i in (0, 2, 4)],
                         [3309, 3310, 3312])
        self.assertEqual([r["tested_budget_passes"] for r in rows], [True, False, False, False, False])


if __name__ == "__main__":
    unittest.main()
