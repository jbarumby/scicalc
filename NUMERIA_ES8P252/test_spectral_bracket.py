"""Adversarial and exact-arithmetic tests for the P252 spectral replays."""

from fractions import Fraction
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P252.odd_pole_audit import certify_odd_pole
from NUMERIA_ES8P252.spectral_bracket import AGGREGATE, bracket_smallest, check_shifted_center, run_audit


def dyadic(pair):
    mantissa, exponent = map(int, pair)
    return Fraction(mantissa) * Fraction(2)**exponent


class SpectralBracketTests(unittest.TestCase):
    def test_exact_known_minimum_with_integration_error(self):
        with ctx.workprec(128):
            center = arb_mat([[3, 1], [1, 3]])
            trial = arb_mat([[1], [-1]])
            error = Fraction(1, 32)
            result = bracket_smallest(center, trial, error)
        self.assertEqual(result["shifted_LDL_positive_pivots"], 2)
        self.assertLess(dyadic(result["lambda_min_strict_lower_dyadic"]), 2 - error)
        self.assertEqual(dyadic(result["lambda_min_upper_dyadic"]), 2 + error)
        self.assertEqual(dyadic(result["quadrature_and_analytic_norm_error_upper_dyadic"]), error)

    def test_shift_rejects_eigenvalue_boundary_and_interval_perturbation(self):
        with ctx.workprec(128):
            center = arb_mat([[3, 1], [1, 3]])
            self.assertEqual(len(check_shifted_center(center, arb(1), arb(0))), 2)
            with self.assertRaisesRegex(ArithmeticError, "LDL"):
                check_shifted_center(center, arb(2), arb(0))
            with self.assertRaisesRegex(ArithmeticError, "LDL"):
                check_shifted_center(center, arb(1), arb(2))

    def test_inexact_or_uninformative_trial_fails_closed(self):
        with ctx.workprec(128):
            center = arb_mat([[3, 1], [1, 3]])
            with self.assertRaisesRegex(ValueError, "exact dyadic"):
                bracket_smallest(center, arb_mat([[arb(1, "1e-20")], [-1]]), Fraction(0))
            with self.assertRaisesRegex(ArithmeticError, "attempted shifts"):
                bracket_smallest(center, arb_mat([[1], [1]]), Fraction(0))

    def test_replays_reject_aggregate_digest_mismatch(self):
        fake_input = (None, None), [], "0" * 64
        with patch("NUMERIA_ES8P252.spectral_bracket.load_aggregate", return_value=fake_input):
            with self.assertRaisesRegex(ArithmeticError, "SHA-256"):
                run_audit(AGGREGATE)
        with patch("NUMERIA_ES8P252.odd_pole_audit.load_aggregate", return_value=fake_input):
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                certify_odd_pole(AGGREGATE)

    def test_archived_exact_bounds_agree_with_the_odd_pole_check(self):
        data = Path(__file__).resolve().parent / "data"
        brackets = json.loads((data / "P252_SPECTRAL_BRACKET.json").read_text(encoding="utf-8"))
        pole = json.loads((data / "P252_ODD_POLE_AUDIT.json").read_text(encoding="utf-8"))
        self.assertEqual(brackets["input_GL60_aggregate_sha256"], pole["input_GL60_aggregate_sha256"])
        for sector in brackets["sectors"].values():
            lower = dyadic(sector["lambda_min_strict_lower_dyadic"])
            upper = dyadic(sector["lambda_min_upper_dyadic"])
            width = dyadic(sector["bracket_width_upper_dyadic"])
            self.assertEqual(sector["shifted_LDL_positive_pivots"], 259)
            self.assertLess(0, lower)
            self.assertLessEqual(upper - lower, width)
            self.assertLess(width, Fraction("5.562e-47"))
        odd = brackets["sectors"]["odd"]
        pole_margin = Fraction(pole["finite_head_margin_strict_lower_exact"])
        self.assertLess(pole_margin, dyadic(odd["lambda_min_strict_lower_dyadic"]))
        self.assertEqual(pole["shifted_pole_free_positive_pivots"], 259)
        scalar_lower, scalar_upper = map(dyadic, pole["one_minus_twice_pole_inverse_pole_endpoints_dyadic"])
        self.assertGreater(scalar_lower, Fraction("7.8412e-45"))
        self.assertLessEqual(scalar_lower, scalar_upper)


if __name__ == "__main__":
    unittest.main()
