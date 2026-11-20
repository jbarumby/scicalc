"""Independent exact arithmetic checks for the finite Gate C support extension."""

from fractions import Fraction
import json
from math import factorial
from pathlib import Path
import unittest

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P251.gate_c.bridge_audit import check_model_floor, local_extension
from NUMERIA_ES8P251.gate_c.spectral_diagnostic import _rayleigh


SYMBOL_FLOOR = Fraction(1794, 10**7)
BASE_HALF_WIDTH = Fraction(161, 200)


def log_upper(integer):
    """Bound log(integer) from above by the positive artanh sequence."""
    ratio = Fraction(integer - 1, integer + 1)
    terms = 16
    series = sum(ratio ** (2 * k + 1) / (2 * k + 1) for k in range(terms))
    first_omitted = ratio ** (2 * terms + 1) / (2 * terms + 1)
    return 2 * (series + first_omitted / (1 - ratio**2))


class LocalExtensionTests(unittest.TestCase):
    def test_exact_right_endpoint_and_remaining_floor(self):
        result = local_extension(SYMBOL_FLOOR)
        increment = Fraction(5, 10**22)
        self.assertEqual(Fraction(result["half_width_right_exact"]), BASE_HALF_WIDTH + increment)
        loss = Fraction(result["common_band_and_pole_loss_upper_exact"])
        remaining = Fraction(result["model_floor_lower_exact"])
        self.assertEqual(loss, Fraction(107, 10**20))
        self.assertEqual(remaining, Fraction(2743, 5 * 10**21))
        self.assertEqual(remaining + loss, Fraction(8093, 5 * 10**21))
        band_floor = Fraction(result["new_band_symbol_lower_exact"])
        self.assertEqual(band_floor, SYMBOL_FLOOR - Fraction(29, 2) * 309 * increment / BASE_HALF_WIDTH)
        self.assertGreater(band_floor, 0)
        derivative = Fraction(result["common_band_and_pole_derivative_upper_exact"])
        self.assertEqual(derivative, Fraction(344364, 161))
        self.assertLess(derivative, loss / increment)

    def test_right_endpoint_stays_within_same_prime_power_segment(self):
        result = local_extension(SYMBOL_FLOOR)
        doubled = 2 * Fraction(result["half_width_right_exact"])
        degree = 18
        partial = sum(doubled**n / factorial(n) for n in range(degree + 1))
        first_omitted = doubled ** (degree + 1) / factorial(degree + 1)
        # The ratios after the first omitted term decrease from 2L/(n+2).
        self.assertGreater(partial, 5)
        self.assertLess(partial + first_omitted / (1 - doubled / (degree + 2)), 6)
        self.assertEqual(result["prime_powers"], [2, 3, 4, 5])
        slope_bound = Fraction(result["prime_slope_upper_exact"])
        self.assertEqual(slope_bound, Fraction(256627, 46750))
        independently_bounded = Fraction(0)
        for power, prime, root in ((2, 2, 1414), (3, 3, 1732), (4, 2, 2000), (5, 5, 2236)):
            root_lower = Fraction(root, 1000)
            self.assertLessEqual(root_lower**2, power)
            independently_bounded += 2 * log_upper(prime) * log_upper(power) / root_lower
        self.assertLess(independently_bounded, slope_bound)

    def test_using_no_increment_preserves_the_fixed_window_floor(self):
        result = local_extension(SYMBOL_FLOOR, half_width_increment=Fraction(0))
        self.assertEqual(Fraction(result["half_width_right_exact"]), BASE_HALF_WIDTH)
        self.assertEqual(Fraction(result["common_band_and_pole_loss_upper_exact"]), 0)
        self.assertEqual(Fraction(result["model_floor_lower_exact"]), Fraction(8093, 5 * 10**21))

    def test_negative_width_or_uncertified_bounds_fail_closed(self):
        with self.assertRaises(ValueError):
            local_extension(SYMBOL_FLOOR, half_width_increment=-Fraction(1, 10**22))
        with self.assertRaisesRegex(ArithmeticError, "frequency band is inconclusive"):
            local_extension(Fraction(0))
        with self.assertRaisesRegex(ArithmeticError, "coercivity"):
            local_extension(SYMBOL_FLOOR, half_width_increment=Fraction(1, 10**20))
        with self.assertRaisesRegex(ArithmeticError, "prime-power window"):
            local_extension(SYMBOL_FLOOR, half_width_increment=1)

    def test_gate_b_must_supply_both_parities_and_strict_reference_margin(self):
        certificate = Path(__file__).resolve().parents[1] / "gate_b/data/P251_GATE_B_CERT.json"
        report = json.loads(certificate.read_text(encoding="utf-8"))
        baseline = Fraction(8093, 5 * 10**21)
        self.assertGreater(check_model_floor(report), baseline)
        report["sectors"]["even"]["certified_positive"] = False
        with self.assertRaisesRegex(ArithmeticError, "both Gate B sectors"):
            check_model_floor(report)
        report["sectors"]["even"]["certified_positive"] = True
        report["schur_transfer"]["model_floor_lower_exact"] = str(baseline)
        with self.assertRaisesRegex(ArithmeticError, "Gate B margin"):
            check_model_floor(report)


class SpectralDiagnosticTests(unittest.TestCase):
    def test_rayleigh_encloses_endpoint_quotients_and_requires_exact_nonzero_vector(self):
        with ctx.workprec(128):
            radius = arb(1) / 8
            matrix = arb_mat([[arb(2, radius), 0], [0, arb(5, radius)]])
            value = _rayleigh(matrix, arb_mat([[1], [2]]))
            center = arb(22) / 5
            for endpoint in (center - radius, center + radius):
                self.assertTrue(value.contains(endpoint))
            for invalid in (arb_mat([[0], [0]]), arb_mat([[arb(1, radius)], [2]])):
                with self.assertRaises(ValueError):
                    _rayleigh(matrix, invalid)


if __name__ == "__main__":
    unittest.main()
