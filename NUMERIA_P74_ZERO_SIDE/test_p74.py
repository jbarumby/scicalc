"""Independent checks for the zero-side interval certificate."""

import copy
from fractions import Fraction
import unittest

from flint import acb, arb, ctx

from NUMERIA_ES8P259.suzuki import SuzukiFormula, prime_powers

from .certify import certify, secant_lower, verify_certificate
from .tails import HEIGHT, TailBounds
from .zeros import verified_ordinates


def rational_ball(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def dyadic_ball(encoded):
    mantissa, exponent = map(int, encoded)
    return rational_ball(mantissa * Fraction(2)**exponent)


def cosine_sum(ordinates, t):
    # Reconstruct the sum using cosines in place of the implementation's sine squares.
    return sum((2 * (1 - (gamma * t).cos()) / gamma**2 for gamma in ordinates), arb(0))


def hardy_sign_at(t):
    term = acb(arb(1) / 2, t).zeta()
    theta = acb(arb(1) / 4, t / 2).lgamma().imag - t * arb.pi().log() / 2
    return term.real * theta.cos() - term.imag * theta.sin()


class VerifiedOrdinateTests(unittest.TestCase):
    def test_order_total_count_and_independent_hardy_signs(self):
        with ctx.workprec(192):
            ordinates, metadata = verified_ordinates(16)
            self.assertEqual(len(ordinates), 16)
            self.assertEqual(metadata["zero_count"], 16)
            self.assertEqual(metadata["count_at_height"], 16)
            self.assertEqual(metadata["hardy_endpoint_sign_verified_count"], 16)
            self.assertGreater(dyadic_ball(metadata["hardy_endpoint_sign_minimum_dyadic"]), arb(0))
            cutoff = dyadic_ball(metadata["count_height_dyadic"])
            self.assertGreater(cutoff, ordinates[-1])
            self.assertEqual(cutoff.zeta_nzeros(), 16)
            self.assertTrue(arb(14) < ordinates[0] < arb(15))
            self.assertTrue(arb(67) < ordinates[-1] < arb(68))
            for prior, following in zip(ordinates, ordinates[1:]):
                self.assertLess(prior, following)
            for gamma in (ordinates[0], ordinates[-1]):
                lo, hi = hardy_sign_at(gamma.lower()), hardy_sign_at(gamma.upper())
                self.assertTrue((lo < 0 < hi) or (hi < 0 < lo))
            self.assertRegex(metadata["ordered_ordinate_bounds_sha256"], r"^[0-9a-f]{64}$")

    def test_precision_overlap(self):
        with ctx.workprec(192):
            first, _ = verified_ordinates(10)
        with ctx.workprec(256):
            second, _ = verified_ordinates(10)
            self.assertTrue(all(a.overlaps(b) for a, b in zip(first, second, strict=True)))


class TailTests(unittest.TestCase):
    def test_coarse_and_graded_monotonicity(self):
        with ctx.workprec(192):
            self.assertEqual(HEIGHT, 3_000_000_000_000)
            coarse, graded = TailBounds("coarse"), TailBounds("graded")
            self.assertTrue(coarse.bound(arb(0)).overlaps(coarse.coarse(arb(0))))
            self.assertLess(coarse.bound(arb(0)), arb("2.96e-12"))
            previous_coarse = previous_graded = arb(0)
            for t in (arb(18), arb(47), arb(55), arb(60)):
                with self.subTest(t=t):
                    coarse_cap, graded_cap = coarse.bound(t), graded.bound(t)
                    self.assertGreater(coarse_cap, previous_coarse)
                    self.assertGreater(graded_cap, previous_graded)
                    self.assertLess(graded_cap, coarse_cap)
                    previous_coarse, previous_graded = coarse_cap, graded_cap

    def test_quartet_loss_covers_sum_with_both_reflections(self):
        with ctx.workprec(192):
            t, gamma, offset = arb(50), arb(HEIGHT + 1), arb(3) / 8
            z = acb(gamma, offset)
            quartet = (z, -z, z.conjugate(), -z.conjugate())
            value = sum(((1 - (root * t).cos()) / root**2 for root in quartet), acb(0))
            self.assertTrue(value.imag.contains(0))
            loss = 4 * (offset * t).cosh() / gamma**2
            self.assertGreater(value.real, -loss)

    def test_normalization_against_existing_prime_formula(self):
        with ctx.workprec(192):
            t = arb(4)
            ordinates, _ = verified_ordinates(32)
            partial_value = cosine_sum(ordinates, t)
            prime_value, _ = SuzukiFormula().values(t)
            self.assertTrue(arb(54) < t.exp() < arb(55))
            for q, p in prime_powers(54):
                prime_value -= arb(p).log() * (t - arb(q).log()) / arb(q).sqrt()
            self.assertGreater(prime_value, partial_value - TailBounds("coarse").bound(t))
            self.assertLess(prime_value - partial_value, arb(1) / 20)


class SecantTests(unittest.TestCase):
    def test_sharp_curvature_loss_at_the_interior(self):
        # On [-1,1], the endpoints of x**2 are both one and the center is zero.
        self.assertEqual(secant_lower(arb(1), arb(1), arb(2), 1), arb(0))

    def test_asymmetric_parabola_switches_between_interior_and_endpoint(self):
        self.assertEqual(secant_lower(arb(3), arb(1), arb(2), 1), arb(3) / 4)
        self.assertEqual(secant_lower(arb(1), arb(7), arb(2), 1), arb(1))

    def test_exported_bound_below_independent_evaluation_on_its_cell(self):
        with ctx.workprec(192):
            record = certify("47.47", "47.48", 16, 192, "graded")
            self.assertEqual(record["status"], "certified_positive")
            left, right = map(Fraction, record["minimum"]["cell"])
            lower = dyadic_ball(record["minimum"]["lower_dyadic"])
            ordinates, _ = verified_ordinates(16)
            tails = TailBounds("graded")
            self.assertGreater(lower, arb(0))
            for weight in (Fraction(0), Fraction(1, 3), Fraction(1, 2), Fraction(1)):
                with self.subTest(weight=weight):
                    t = rational_ball(left + weight * (right - left))
                    self.assertGreater(cosine_sum(ordinates, t) - tails.bound(t), lower)


class CertificateTests(unittest.TestCase):
    def test_exact_rational_endpoint_and_replay(self):
        first = certify("47.47", "47.48", 16, 192, "graded")
        same = certify(Fraction(4747, 100), Fraction(1187, 25), 16, 192, "graded")
        self.assertEqual(first, same)
        self.assertEqual(first["end_exact"], "1187/25")
        self.assertEqual(first["certified_through_exact"], "1187/25")
        self.assertEqual(first["tail_profile"], "graded")
        self.assertGreater(first["cell_count"], 0)
        self.assertIsNone(first["first_unresolved_cell"])
        verify_certificate(first)

    def test_inconclusive_interval_cannot_export_positive_coverage(self):
        result = certify("50", "50.015625", 16, 192, "coarse")
        self.assertEqual(result["status"], "inconclusive")
        self.assertEqual(result["cell_count"], 0)
        self.assertIsNone(result["minimum"])
        self.assertIsNotNone(result["first_unresolved_cell"])
        self.assertNotEqual(result["certified_through_exact"], result["end_exact"])
        verify_certificate(result)

    def test_rigorous_uniform_and_adaptive_methods(self):
        for method in ("uniform", "adaptive"):
            with self.subTest(method=method):
                result = certify(18, "18.03125", 16, 192, "coarse", method=method)
                self.assertEqual(result["method"], method)
                self.assertEqual(result["status"], "certified_positive")
                self.assertEqual(result["certified_through_exact"], result["end_exact"])

    def test_checkpoint_replays_cumulative_coverage(self):
        result = certify(18, "18.03", 16, 192, "coarse", checkpoints=("18.01", "18.02"))
        self.assertEqual(result["status"], "certified_positive")
        rows = result["checkpoints"]
        self.assertEqual([row["target_exact"] for row in rows], ["1801/100", "901/50", "1803/100"])
        self.assertEqual([row["cells_through_target"] for row in rows], [1, 2, 3])
        with ctx.workprec(192):
            bounds = [dyadic_ball(row["minimum_through_target_dyadic"]) for row in rows]
            self.assertTrue(bounds[0] >= bounds[1] >= bounds[2] > 0)
        verify_certificate(result)

    def test_reject_invalid_and_nonexact_inputs(self):
        for start, end, count, precision in (
            (18, 18, 16, 192),
            (19, 18, 16, 192),
            (-1, 18, 16, 192),
            (18.0, 19, 16, 192),
            (18, 19, 0, 192),
            (18, 19, True, 192),
            (18, 19, 16, 0),
        ):
            with self.subTest(start=start, end=end, count=count, precision=precision):
                with self.assertRaises(ValueError):
                    certify(start, end, count, precision)

    def test_reject_tampered_certificates_and_type_aliases(self):
        original = certify(18, "18.03125", 16, 192, "coarse")
        for field, value in (
            ("schema", "bogus"),
            ("irrelevant_field", True),
            ("status", "inconclusive"),
            ("cell_count", float(original["cell_count"])),
            ("certified_through_exact", "18"),
        ):
            with self.subTest(field=field):
                damaged = copy.deepcopy(original)
                damaged[field] = value
                with self.assertRaises(ValueError):
                    verify_certificate(damaged)
        damaged = copy.deepcopy(original)
        damaged["minimum"]["lower_dyadic"][0] = "1"
        with self.assertRaises(ValueError):
            verify_certificate(damaged)


if __name__ == "__main__":
    unittest.main()
