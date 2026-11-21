"""Regression checks for scalar Suzuki and the banded Weil auxiliary."""

from fractions import Fraction
import unittest

from flint import arb, ctx

from NUMERIA_ES8P258.concentration import analytic_concentration, compare_concentration, concentration_chunk
from .banded_weil import certify_analytic_constants, chebyshev_tail_bounds, effective_head_step
from .suzuki import SuzukiFormula, certify, prime_powers


class SuzukiTest(unittest.TestCase):
    def test_sieve_includes_powers_once(self):
        self.assertEqual(prime_powers(16), [(2, 2), (3, 3), (4, 2), (5, 5), (7, 7),
                                          (8, 2), (9, 3), (11, 11), (13, 13), (16, 2)])

    def test_initial_interval(self):
        with ctx.workprec(192):
            self.assertGreater(SuzukiFormula().initial_margin(), Fraction(1, 4))

    def test_derivative_against_closed_form(self):
        with ctx.workprec(192):
            formula = SuzukiFormula()
            for t in (arb(1) / 64, arb(2).log(), arb(7), arb(1000).log()):
                _, slope = formula.values(t)
                r = (-t / 2).exp()
                independent = 2 * (1 / r - r) + formula.gamma_coefficient + r.atanh() + r.atan()
                self.assertTrue(slope.overlaps(independent))

    def test_bound_at_prime_power_contact(self):
        with ctx.workprec(192):
            formula = SuzukiFormula()
            t = arb(7).log()
            raw, _ = formula.values(t)
            before = raw
            for q, p in prime_powers(6):
                before -= arb(p).log() / arb(q).sqrt() * (t - arb(q).log())
            after = before - arb(7).log() / arb(7).sqrt() * (t - arb(7).log())
            self.assertTrue(before.overlaps(after))
            self.assertTrue(before > 0)

    def test_full_small_continuous_certificates(self):
        self.assertEqual(certify(2)["positive_for_all_t_in"], "(0,log(2)]")
        report = certify(128)
        self.assertEqual(report["positive_for_all_t_in"], "(0,log(128)]")
        self.assertGreater(report["certified_cells_after_initial_interval"], 25)
        self.assertEqual(report["all_affine_prefixes_verified_through_limit"], 128)
        self.assertFalse(report["unbounded_range_established"])


class BandedWeilTest(unittest.TestCase):
    def test_strip_and_complement_constants(self):
        self.assertEqual(certify_analytic_constants()["strip_loss_integer_upper"], 1839)
        self.assertEqual(effective_head_step(Fraction(1, 4), Fraction(1, 8)), Fraction(1, 12))
        with self.assertRaises(ValueError):
            effective_head_step(Fraction(1, 4), 0)

    def test_advanced_band_has_a_contiguous_adaptive_proof(self):
        cert = certify_analytic_constants(590)["high_band_adaptive_certificate"]
        self.assertEqual(cert["adaptive_range"], [590, 850])
        self.assertEqual(cert["adaptive_cover"][0][0], "590")
        self.assertEqual(cert["adaptive_cover"][-1][1], "850")
        for left, right in zip(cert["adaptive_cover"], cert["adaptive_cover"][1:]):
            self.assertEqual(left[1], right[0])

    def test_direct_chebyshev_bounds_for_590(self):
        for parity in (0, 1):
            coupling_squared, deviation = chebyshev_tail_bounds(590, 377, parity)
            self.assertTrue(0 < coupling_squared < Fraction("1e-65"))
            self.assertTrue(0 < deviation < Fraction("1e-65"))
        with self.assertRaises(ValueError):
            chebyshev_tail_bounds(850, 377, 0)

    def test_G850_endpoint_formula_against_a_direct_GL60_panel(self):
        args = {"cutoff": 850, "modes": 4}
        q = concentration_chunk(1699, 1700, **args)
        h = analytic_concentration(1699, 1700, **args)
        report = compare_concentration(q, h, **args)
        self.assertEqual([r["endpoint_enclosures_contained"] for r in report.values()], [10, 10])


if __name__ == "__main__":
    unittest.main()
