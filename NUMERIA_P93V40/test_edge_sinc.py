"""Independent integration tests for the sinc endpoint representer."""

import unittest
from fractions import Fraction

from flint import acb, arb, ctx

from .edge_sinc import endpoint_coefficients, endpoint_norm_squares, endpoint_tail_norm_squares


class SincEdgeTests(unittest.TestCase):
    def test_legendre_coefficients_against_direct_acb_quadrature(self):
        width, omega = Fraction(3, 2), 3
        with ctx.workprec(256):
            a = arb(3) / 2
            actual = endpoint_coefficients(5, width=width, omega=omega)
            for n in range(10):
                direct = acb.integral(
                    lambda x, _: ((a-x)*omega).sinc() * (x/a).legendre_p(n),
                    -a, a, rel_tol=arb("1e-65"), abs_tol=arb("1e-65"),
                )
                direct *= omega / arb.pi() * ((2*n+1) / (2*a)).sqrt()
                with self.subTest(degree=n):
                    self.assertTrue(direct.imag.contains(0))
                    self.assertTrue(actual[n % 2][n // 2].overlaps(direct.real))
                    self.assertLess(actual[n % 2][n // 2].rad(), arb("1e-68"))

    def test_parity_norms_against_direct_acb_quadrature(self):
        width, omega = Fraction(3, 2), 3
        with ctx.workprec(256):
            a = arb(3) / 2
            norms = endpoint_norm_squares(width=width, omega=omega)
            coefficients = endpoint_coefficients(5, width=width, omega=omega)
            tails = endpoint_tail_norm_squares(coefficients, width=width, omega=omega)
            for parity in (0, 1):
                direct = acb.integral(
                    lambda x, _: (((a-x)*omega).sinc() + (-1)**parity*((a+x)*omega).sinc())**2,
                    -a, a, rel_tol=arb("1e-65"), abs_tol=arb("1e-65"),
                ) * (omega / (2 * arb.pi()))**2
                with self.subTest(parity=parity):
                    self.assertTrue(direct.imag.contains(0))
                    self.assertTrue(norms[parity].overlaps(direct.real))
                    self.assertGreater(tails[parity], 0)
                    self.assertLess(tails[parity], arb("1e-7"))

    def test_stable_at_physical_cutoff(self):
        with ctx.workprec(256):
            short = endpoint_coefficients(4)
            long = endpoint_coefficients(200)
            full = endpoint_norm_squares()
            tail = endpoint_tail_norm_squares(long)
            for p in (0, 1):
                self.assertTrue(all(a.overlaps(b) for a, b in zip(short[p], long[p])))
                self.assertLess(max(value.rad() for value in long[p]), arb(2)**-230)
                self.assertGreater(tail[p], 0)
                self.assertLess(tail[p], full[p])

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            endpoint_coefficients(0)
        with self.assertRaises(ValueError):
            endpoint_coefficients(2, omega=-1)
        with self.assertRaises(ValueError):
            endpoint_norm_squares(width=0)


if __name__ == "__main__":
    unittest.main()
