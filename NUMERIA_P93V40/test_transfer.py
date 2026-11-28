"""Tests for the conditional Schur and prime-16 strip estimates."""

from fractions import Fraction
import unittest

from .transfer import certify_transfer


def _exact_dyadic(encoded):
    mantissa, exponent = encoded
    return int(mantissa) * Fraction(2)**exponent


class TransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = certify_transfer()

    def test_width_crosses_log16_by_independent_rational_series(self):
        # log(2) = 2 atanh(1/3). Bound the series remainder geometrically.
        terms = 92
        log2_lower = sum(
            Fraction(2, (2 * k + 1) * 3**(2 * k + 1)) for k in range(terms)
        )
        log2_upper = log2_lower + Fraction(9, 4 * (2 * terms + 1) * 3**(2 * terms + 1))
        width = Fraction(self.report["width_exact"])
        independent_lower = 2 * width - 4 * log2_upper
        independent_upper = 2 * width - 4 * log2_lower
        reported_lower = _exact_dyadic(self.report["delta_lower_dyadic"])
        reported_upper = _exact_dyadic(self.report["delta_upper_dyadic"])
        self.assertGreater(independent_lower, 0)
        self.assertLess(independent_upper, Fraction(3, 10**75))
        self.assertGreater(reported_lower, 0)
        self.assertLessEqual(reported_lower, independent_upper)
        self.assertGreaterEqual(reported_upper, independent_lower)

    def test_transfer_spends_less_than_the_head_margin(self):
        h = Fraction(1, 10**60)
        delta_upper = _exact_dyadic(self.report["delta_upper_dyadic"])
        loss = _exact_dyadic(self.report["new_prime_loss_upper_dyadic"])
        retained = _exact_dyadic(self.report["conditional_crossing_margin_lower_dyadic"])
        self.assertIs(self.report["conditional_on_complete_heads_and_spectral_majorant"], True)
        self.assertEqual(Fraction(self.report["head_target_exact"]), h)
        self.assertGreater(_exact_dyadic(self.report["far_symbol_lower_dyadic"]), 3)
        self.assertLess(_exact_dyadic(self.report["low_abs_symbol_upper_dyadic"]), 18)
        self.assertGreaterEqual(loss, Fraction(4_000_000, 6) * delta_upper)
        self.assertGreater(retained, Fraction(89, 100) * h)
        self.assertLessEqual(retained + loss, Fraction(9, 10) * h)

    def test_separated_strips_improve_the_uniform_budget(self):
        original_loss = _exact_dyadic(self.report["new_prime_loss_upper_dyadic"])
        sharper = certify_transfer("8e-70", separated_strips=True)
        loss = _exact_dyadic(sharper["new_prime_loss_upper_dyadic"])
        retained = _exact_dyadic(sharper["conditional_crossing_margin_lower_dyadic"])
        self.assertEqual(sharper["strip_cutoff"], 1_600_000)
        self.assertGreater(_exact_dyadic(sharper["strip_bound_formula_lower_dyadic"]),
                           Fraction("1.65e-70"))
        self.assertLess(loss, original_loss / 4)
        self.assertGreater(retained, Fraction("5.54e-70"))
        with self.assertRaisesRegex(ArithmeticError, "margin is exhausted"):
            certify_transfer("1e-70", separated_strips=True)

    def test_ellipse_l2_bound_brings_the_five_schur_gates_under_3310(self):
        result = certify_transfer("8e-70", separated_strips=True, chebyshev_tail=True,
                                  retention="1/4")
        self.assertEqual(Fraction(result["retention_exact"]), Fraction(1, 4))
        self.assertGreater(_exact_dyadic(result["conditional_crossing_margin_lower_dyadic"]),
                           Fraction("3.47e-71"))
        self.assertEqual([row["symbol_cap_18"]["first_sufficient_dimension"]
                          for row in result["theta_results"]], list(range(3305, 3310)))

    def test_all_five_schur_dimensions_at_fixed_rho(self):
        self.assertEqual(Fraction(self.report["rho_exact"]), Fraction(33, 25))
        rows = self.report["theta_results"]
        self.assertEqual([Fraction(row["theta"]) for row in rows],
                         [Fraction(i, 4) for i in range(5)])
        for label, sizes in (
            ("symbol_cap_50", [3309, 3310, 3310, 3311, 3312]),
            ("symbol_cap_18", [3307, 3308, 3309, 3310, 3311]),
        ):
            for row, size in zip(rows, sizes, strict=True):
                with self.subTest(theta=row["theta"], cap=label):
                    gate = row[label]
                    self.assertEqual(gate["first_sufficient_dimension"], size)
                    self.assertEqual(gate["passes_at_3309"], size <= 3309)
                    self.assertGreater(_exact_dyadic(gate["g_even_upper_dyadic"]), 0)
                    self.assertGreater(_exact_dyadic(gate["tail_reserve_lower_dyadic"]), 0)
                    self.assertGreater(_exact_dyadic(gate["schur_slack_lower_dyadic"]), 0)


if __name__ == "__main__":
    unittest.main()
