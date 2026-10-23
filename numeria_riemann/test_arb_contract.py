"""Reproducible python-flint 0.9.0 Arb/Acb contract tests.

The production Weil audit uses real Arb range-sum quadrature. These tests
exercise Acb's separate complex-segment integration API and the removable
sinc singularity without changing the certificate algorithm.
"""
import math
import unittest

import flint
from flint import acb, arb, ctx

from numeria_riemann.weil_interval_audit import _sinc


def _sinc_acb_with_removable_extension(z, analytic):
    """Enclose sin(z)/z on complex balls, including balls containing zero.

    On |z| <= 1 use the entire Taylor series through z**64/65! and a
    rectangular error enclosure of +/-1e-90 in each component. The true
    absolute tail is < 1e-90 for |z| <= 1. Outside that disk, a ball
    containing zero is conservatively returned as non-finite.
    """
    if z.contains(0):
        if not z.abs_upper() < arb(1):
            return acb(arb("nan"), arb("nan"))
        total = acb(1)
        term = acb(1)
        for n in range(1, 33):
            term = term * z * z / ((2 * n) * (2 * n + 1))
            total += term if n % 2 == 0 else -term
        error = arb("0 +/- 1e-90")
        return total + acb(error, error)
    return z.sin() / z


class ArbAcbContractTests(unittest.TestCase):
    def setUp(self):
        self.old_precision = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self.old_precision

    def test_fixed_python_flint_version(self):
        self.assertEqual(flint.__version__, "0.9.0")

    def test_real_sinc_at_exact_removable_singularity(self):
        value = _sinc(arb(0))
        self.assertTrue(value.is_finite())
        self.assertTrue(value.contains(1))

    def test_real_sinc_ball_containing_zero_is_enclosed(self):
        value = _sinc(arb(0, "0.01"))
        self.assertTrue(value.is_finite())
        self.assertTrue(value.contains(1))

    def test_real_sinc_rejects_too_wide_ball_crossing_zero(self):
        with self.assertRaises(ValueError):
            _sinc(arb(0, "0.2"))

    def test_acb_integral_contract_on_complex_segment(self):
        # Integral of z dz from 0 to 1+i is ((1+i)^2)/2 = i.
        value = acb.integral(lambda z, analytic: z, acb(0), acb(1, 1))
        self.assertTrue(value.is_finite())
        self.assertTrue(value.real.contains(0))
        self.assertTrue(value.imag.contains(1))

    def test_acb_integral_of_sinc_crossing_zero(self):
        # The exact integral is 2*Si(1/20); enclose it independently by
        # integrating the even Taylor series term by term.
        value = acb.integral(
            _sinc_acb_with_removable_extension,
            acb("-1/20"),
            acb("1/20"),
        )
        self.assertTrue(value.is_finite())
        self.assertTrue(value.imag.contains(0))

        a = arb("1/20")
        expected = arb(0)
        for n in range(33):
            term = 2 * a ** (2 * n + 1) / (
                (2 * n + 1) * math.factorial(2 * n + 1)
            )
            expected += term if n % 2 == 0 else -term
        expected += arb("0 +/- 1e-90")
        self.assertTrue(value.real.overlaps(expected))


if __name__ == "__main__":
    unittest.main()
