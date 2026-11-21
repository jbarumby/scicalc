"""Independent checks for the segmented sieve and Suzuki prefix certificate."""

import copy
from fractions import Fraction
from math import isqrt
import unittest
from unittest.mock import patch

from flint import arb, ctx

from NUMERIA_ES8P259.suzuki import SuzukiFormula, certify as certify_p259, prime_powers as powers_p259

from .certify import certify, quadratic_gap, tangent_gap, verify_certificate
from .prime_powers import iter_prime_powers


def trial_powers(limit):
    """Identify the prime powers by removing the smallest divisor of each integer."""
    for q in range(2, limit + 1):
        p = next((d for d in range(2, isqrt(q) + 1) if q % d == 0), q)
        remaining = q
        while remaining % p == 0:
            remaining //= p
        if remaining == 1:
            yield q, p


def log_bounds(x):
    """Use the rational series for 2*atanh((x-1)/(x+1)) with a geometric tail."""
    y = (x - 1) / (x + 1)
    n = 48
    total = sum((2 * y ** (2 * k + 1) / (2 * k + 1) for k in range(n)), Fraction(0))
    tail = 2 * y ** (2 * n + 1) / ((2 * n + 1) * (1 - y * y))
    return total, total + tail


def sqrt_bounds(q):
    scale = 2**64
    lower = isqrt(q * scale * scale)
    return Fraction(lower, scale), Fraction(lower + 1, scale)


def dyadic_fraction(pair):
    mantissa, exponent = map(int, pair)
    return mantissa * Fraction(2) ** exponent


class PrimePowerTests(unittest.TestCase):
    def test_trial_division_across_segment_boundaries(self):
        for limit in (2, 3, 7, 8, 9, 25, 26, 27, 64, 121, 512):
            expected = list(trial_powers(limit))
            for block_size in (1, 2, 3, 8, 9, 31, 64, 1024):
                with self.subTest(limit=limit, block_size=block_size):
                    self.assertEqual(list(iter_prime_powers(limit, block_size)), expected)

    def test_match_original_p259_sieve(self):
        powers = list(iter_prime_powers(10000, 509))
        self.assertEqual(len(powers), 1280)
        self.assertEqual(powers, powers_p259(10000))

    def test_no_powers_beyond_inclusive_limit(self):
        self.assertEqual(list(iter_prime_powers(81, 80))[-1], (81, 3))
        self.assertEqual(list(iter_prime_powers(80, 80))[-1], (79, 79))


class GapTests(unittest.TestCase):
    def test_exact_tangent_and_contact(self):
        with ctx.workprec(192):
            self.assertEqual(tangent_gap(arb(2), arb(7)), arb(11))
            for a, b in ((arb(1) / 3, arb(2)), (arb(7), arb(-3))):
                contact = 2 * (a / 2).log()
                at_contact = 4 * (contact / 2).exp() - a * contact + b
                gap = tangent_gap(a, b)
                self.assertTrue(gap.overlaps(at_contact))
                nearby = contact + arb(1) / 16
                self.assertGreater(4 * (nearby / 2).exp() - a * nearby + b, gap)

    def test_quadratic_bound_against_pure_rational_sums(self):
        with ctx.workprec(192):
            for q in (2, 4, 9):
                with self.subTest(q=q):
                    a_lo, a_hi, p_lo, p_hi = 4 * (Fraction(0),)
                    a, b = arb(0), arb(0)
                    for n, p in trial_powers(q):
                        log_lo, log_hi = log_bounds(Fraction(p))
                        root_lo, root_hi = sqrt_bounds(n)
                        w_lo, w_hi = log_lo / root_hi, log_hi / root_lo
                        log_ratio_lo, log_ratio_hi = log_bounds(Fraction(q, n))
                        a_lo += w_lo
                        a_hi += w_hi
                        p_lo += w_lo * log_ratio_lo
                        p_hi += w_hi * log_ratio_hi
                        w = arb(p).log() / arb(n).sqrt()
                        a += w
                        b += w * arb(n).log()
                    s_lo, s_hi = sqrt_bounds(q)
                    d_lo, d_hi = 4 * s_lo - p_hi, 4 * s_hi - p_lo
                    deviation_lo, deviation_hi = a_lo - 2 * s_hi, a_hi - 2 * s_lo
                    self.assertLess(deviation_hi, 0)
                    q_lo = d_lo - deviation_lo**2 / s_lo
                    q_hi = d_hi - deviation_hi**2 / s_hi

                    direct, lower_bound = quadratic_gap(a, b, q)
                    for actual, lo, hi in ((direct, d_lo, d_hi), (lower_bound, q_lo, q_hi)):
                        self.assertLessEqual(dyadic_fraction(actual.lower().man_exp()), hi)
                        self.assertGreaterEqual(dyadic_fraction(actual.upper().man_exp()), lo)
                    self.assertGreater(q_lo, 0)
                    self.assertGreater(tangent_gap(a, b), lower_bound)


class CertificateTests(unittest.TestCase):
    def test_segment_size_independence_and_p259_provenance(self):
        result = certify(128, block_size=1)
        self.assertEqual(result, certify(128, block_size=43))
        self.assertEqual(result["prime_power_count"], len(list(trial_powers(128))))
        self.assertEqual(
            result["prefix_tangent_certificate_sha256"],
            certify_p259(128)["prefix_tangent_certificate_sha256"],
        )
        for field in ("tangent_minimum", "quadratic_minimum"):
            self.assertEqual(result[field]["prime_power"], 2)
            self.assertGreater(dyadic_fraction(result[field]["lower_dyadic"]), 1)
        self.assertIsNone(verify_certificate(result, block_size=19))

    def test_overlap_includes_a_composite_extension_endpoint(self):
        with patch("NUMERIA_SUZUKI_PREFIX.certify.P259_LIMIT", 9):
            original = certify(9)
            extension = certify(10)
        self.assertEqual(extension["p259_overlap_limit"], 9)
        self.assertEqual(extension["p259_overlap_sha256"], original["prefix_tangent_certificate_sha256"])
        self.assertEqual(extension["new_prime_power_count"], 0)
        self.assertEqual(extension["screw_extension"]["active_prefix_count"], 1)

    def test_screw_extension_bound_against_p259_lerch_series(self):
        with patch("NUMERIA_SUZUKI_PREFIX.certify.P259_LIMIT", 64):
            result = certify(128)
        screw = result["screw_extension"]
        self.assertTrue(screw["strictly_positive"])
        self.assertEqual(screw["active_prefix_count"], result["new_prime_power_count"] + 1)
        bound = dyadic_fraction(screw["lower_bound_minimum"]["lower_dyadic"])
        with ctx.workprec(192):
            formula = SuzukiFormula()
            for x in (64, Fraction(129, 2), 81, 127, 128):
                x = Fraction(x)
                t = (arb(x.numerator) / x.denominator).log()
                value, _ = formula.values(t)
                for q, p in trial_powers(int(x)):
                    value -= arb(p).log() / arb(q).sqrt() * (t - arb(q).log())
                self.assertGreater(value, arb(bound.numerator) / bound.denominator)

    def test_reject_tampering_in_values_and_types(self):
        result = certify(2)
        for field, value in (
            ("limit", 3),
            ("prime_power_count", True),
            ("unexpected_field", 0),
            ("prefix_tangent_certificate_sha256", "0" * 64),
            ("quadratic_minimum", {"prime_power": 2, "lower_dyadic": ["7", 0]}),
        ):
            with self.subTest(field=field):
                damaged = copy.deepcopy(result)
                damaged[field] = value
                with self.assertRaises(ValueError):
                    verify_certificate(damaged)
        del result["limit"]
        with self.assertRaises(ValueError):
            verify_certificate(result)

    def test_reject_invalid_precision(self):
        for precision in (0, 63, 192.0):
            with self.subTest(precision=precision), self.assertRaises(ValueError):
                certify(2, precision=precision)


if __name__ == "__main__":
    unittest.main()
