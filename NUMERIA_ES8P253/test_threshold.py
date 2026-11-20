import json
from fractions import Fraction
from pathlib import Path
import unittest

from NUMERIA_ES8P249.legendre_tail import fixed_symbol_half_width, four_prime_model_half_width
from NUMERIA_ES8P251.gate_b.finite_head import _dyadic_fraction
from NUMERIA_ES8P253.prime7_tail import certify_four_prime_floor, certify_tail
from NUMERIA_ES8P253.threshold_audit import AUX_WIDTH, crossing_bounds


DATA = Path(__file__).with_name("data")


class TestPrime7Threshold(unittest.TestCase):
    def test_old_weil_wrapper_does_not_accept_auxiliary_width(self):
        self.assertEqual(four_prime_model_half_width(AUX_WIDTH), AUX_WIDTH)
        with self.assertRaises(ValueError):
            fixed_symbol_half_width(AUX_WIDTH)

    def test_negative_or_outside_auxiliary_data_rejected(self):
        with self.assertRaises(ValueError):
            crossing_bounds(Fraction(0))
        with self.assertRaises(ValueError):
            crossing_bounds(Fraction("1e-29"), half_width=Fraction(19, 20))

    def test_crossing_spends_at_most_half_the_auxiliary_floor(self):
        for floor in (Fraction("1e-27"), Fraction("1e-37")):
            with self.subTest(floor=floor):
                result = crossing_bounds(floor)
                strip = Fraction(result["strip_width_exact"])
                limit = Fraction(result["auxiliary_span_headroom_lower_exact"])
                certified = Fraction(result["Weil_floor_lower_exact"])
                self.assertTrue(0 < strip < limit)
                self.assertEqual(Fraction(result["half_width_increment_above_log7_over_2_exact"]), strip / 2)
                self.assertEqual(certified, floor - 2520 * strip)
                self.assertGreaterEqual(certified, floor / 2)

    def test_signed_tail_has_contiguous_positive_cover_and_negative_witness(self):
        report = certify_tail()
        self.assertEqual(report, json.loads((DATA / "PRIME7_TAIL.json").read_text()))
        end = Fraction(1551)
        for left, right, bound in report["adaptive_cover"]:
            self.assertEqual(Fraction(left), end)
            end = Fraction(right)
            self.assertGreater(end, Fraction(left))
            self.assertGreater(_dyadic_fraction(bound), 0)
        self.assertEqual(end, 2200)
        self.assertGreaterEqual(Fraction(report["rational_negative_witness_frequency"]), 1550)
        self.assertLess(_dyadic_fraction(report["negative_witness_dyadic_endpoints"][1]), 0)

    def test_raised_four_prime_floor(self):
        report = certify_four_prime_floor()
        self.assertEqual(report, json.loads((DATA / "FOUR_PRIME_RAISED_FLOOR.json").read_text()))
        self.assertEqual(report["proposed_beta_exact"], "1/4")
        self.assertEqual(report["certified_range_positive_frequency"], [504, "infinity"])
        self.assertGreater(_dyadic_fraction(report["adaptive_uniform_strict_lower_dyadic"]), 0)
        self.assertGreater(_dyadic_fraction(report["far_tail_uniform_strict_lower_dyadic"]), 0)

    def test_rayleigh_witnesses_for_the_low_cutoffs(self):
        bounds = (
            (309, ("-1.85e-26", "-1.79e-17")),
            (350, ("-4.70e-29", "-8.68e-20")),
            (504, ("-1.12e-29", "-5.30e-23")),
        )
        for cutoff, limits in bounds:
            report = json.loads((DATA / f"T{cutoff}_NEGATIVE_WITNESS.json").read_text())
            self.assertEqual(report["cutoff"], cutoff)
            for sector, upper in zip(("even", "odd"), limits, strict=True):
                with self.subTest(cutoff=cutoff, sector=sector):
                    witness = report["sectors"][sector]
                    self.assertLess(_dyadic_fraction(witness["continuous_Rayleigh_upper_dyadic"]), Fraction(upper))
                    self.assertEqual(witness["vector_nonpadding_coordinates_exact_dyadic"][-1], ["1", 0])


if __name__ == "__main__":
    unittest.main()
