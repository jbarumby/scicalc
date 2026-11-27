"""Independent Arb recurrences and exact-input controls for Binet's budget."""

from fractions import Fraction
from unittest import TestCase

from flint import acb, arb, ctx

from .digamma_binet import PRIME_POWERS, audit, binet_majorants, replay_reserve


def dyadic(point):
    mantissa, exponent = map(int, point)
    return Fraction(mantissa) * Fraction(2)**exponent


class DigammaAuditTests(TestCase):
    def test_binet_and_arb_at_cutoff(self):
        for precision in (192, 256):
            with self.subTest(precision=precision):
                record = audit(precision=precision)
                floor = dyadic(record["uniform_archimedean_binet_floor_lower_dyadic"])
                direct = dyadic(record["uniform_archimedean_digamma_floor_lower_dyadic"])
                loss = dyadic(record["binet_exponential_loss_upper_dyadic"])
                self.assertGreater(floor, Fraction("6.57395558"))
                self.assertGreater(direct, floor)
                self.assertLess(direct - floor, Fraction(1, 10**7))
                self.assertTrue(0 < loss < Fraction(1, 10**3000))

    def test_upper_loss_constants_and_direct_values(self):
        result = audit()
        cap = Fraction(result["refined_upper_loss_constant_exact"])
        self.assertEqual(cap, Fraction("6.1059"))
        self.assertLess(dyadic(result["historical_budget_derived_constant_upper_dyadic"]), 8)
        with ctx.workprec(256):
            for t in (4500, 7000, -4500, -9000, 10**6):
                with self.subTest(t=t):
                    z = arb(t)
                    value = acb(arb(1) / 4, z / 2).digamma().real - arb.pi().log() - 2
                    for q, p in PRIME_POWERS:
                        value -= 2 * arb(p).log() / arb(q).sqrt() * (arb(q).log() * z).cos()
                    self.assertLess(value, abs(z).log() + arb(cap.numerator) / cap.denominator)
        self.assertGreater(dyadic(result["no_deficit_above_beta_lower_dyadic"]), Fraction("0.0338"))

    def test_anchor_and_loss_change_monotonically(self):
        with ctx.workprec(192):
            first, first_loss, first_compact = binet_majorants(4500)
            second, second_loss, second_compact = binet_majorants("4501/1")
            self.assertGreater(second, first)
            self.assertLess(second_loss, first_loss)
            self.assertLess(second_compact, first_compact)

    def test_exact_component_replay(self):
        record = replay_reserve()
        self.assertEqual(record["components"], 1230)
        self.assertEqual(record["multipliers"], {"2": 1110, "14": 120})
        self.assertEqual(Fraction(record["adaptive_spectral_reserve_exact"]),
                         Fraction("2") - Fraction("34164906632832863/18157142016000000"))
        self.assertGreater(Fraction(record["harmonic_spectral_reserve_exact"]), Fraction("0.02"))

    def test_no_inexact_or_out_of_domain_inputs(self):
        for cutoff in (4500.0, True, 0, "-4500", "1/2", "1/0"):
            with self.subTest(cutoff=cutoff):
                with self.assertRaises(ValueError):
                    binet_majorants(cutoff)
        with self.assertRaises(ValueError):
            audit(beta=2.0)
        with self.assertRaises(ValueError):
            audit(precision=64)
