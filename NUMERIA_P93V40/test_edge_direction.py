import unittest

from flint import acb, arb, ctx

from .assembly import WIDTH
from .edge_direction import BAND_EDGES, band_tail_concentration_upper, directional_cross_majorant
from .transfer import ball, pole_audit


class DirectionalCrossMajorantTests(unittest.TestCase):
    def test_band_trace_bound_against_direct_integrals(self):
        with ctx.workprec(256):
            width = ball(WIDTH)
            split, left, right = 3, "0.5", "1.25"
            for parity in (0, 1):
                def density(t, analytic):
                    x = width * t
                    sine, cosine = x.sin(), x.cos()
                    bessel = [sine / x, (sine - x * cosine) / x**2]
                    orders = range(2 * split + parity, 2 * split + parity + 6, 2)
                    for n in range(1, orders[-1]):
                        bessel.append((2 * n + 1) / x * bessel[-1] - bessel[-2])
                    return 2 * width / arb.pi() * sum((2 * n + 1) * bessel[n]**2 for n in orders)

                part_trace = acb.integral(density, acb(arb(left)), acb(arb(right)),
                                          rel_tol=arb("1e-45"))
                bound = band_tail_concentration_upper(split, parity, left, right)
                self.assertTrue(bound.is_exact())
                self.assertTrue(part_trace.imag.contains(0))
                self.assertGreater(part_trace.real, 0)
                self.assertGreater(bound, part_trace.real)
                self.assertLess(bound, 1)

    def test_weighted_endpoints_against_constant_mode_primitive(self):
        with ctx.workprec(320):
            n, width = 3309, ball(WIDTH)
            result = directional_cross_majorant(1, n, 0)[0, 0]

            def primitive(endpoint):
                if endpoint == 0:
                    return arb(0)
                x = width * endpoint
                return 2 / arb.pi() * ((2 * x).si() - x.sin()**2 / x)

            direct = sum(band_tail_concentration_upper(n, 0, l, r) *
                         (primitive(r) - primitive(l)) for l, r in zip(BAND_EDGES, BAND_EDGES[1:]))
            pole, _ = pole_audit(n)
            expected = len(BAND_EDGES) * (18**2 * direct + ball(pole**2))
            self.assertTrue(result.overlaps(expected))
            whole_band = band_tail_concentration_upper(n, 0, BAND_EDGES[0], BAND_EDGES[-1])
            self.assertGreater(result, 0)
            self.assertLess(result, 18**2 * whole_band / 100)

    def test_invalid_prefix_and_band_are_rejected(self):
        with self.assertRaises(ValueError):
            directional_cross_majorant(10, 7, 0)
        with self.assertRaises(ValueError):
            directional_cross_majorant(1, 10, 2)
        with self.assertRaises(ValueError):
            band_tail_concentration_upper(10, 1, -1, 1)


if __name__ == "__main__":
    unittest.main()
