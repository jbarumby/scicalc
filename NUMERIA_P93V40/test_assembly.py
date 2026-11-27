"""Compare the archived-band endpoint formula against direct Arb integrals."""

from fractions import Fraction
import unittest

from flint import acb, arb, ctx

from .assembly import WIDTH, _bessels, band_gram, deficit_endpoints


class BandTests(unittest.TestCase):
    def test_guarded_bessel_encloses_point_and_width(self):
        with ctx.workprec(384):
            x = arb(200, arb(2)**-370)
            values = _bessels(x, 399)
            self.assertLess(max(value.rad() for value in values), arb(2)**-330)
            for point in (x.lower(), x.mid(), x.upper()):
                self.assertTrue(values[0].overlaps(point.sin() / point))
                with ctx.workprec(768):
                    j399 = (arb.pi() / (2*point)).sqrt() * point.bessel_j(arb(799)/2)
                self.assertTrue(values[-1].overlaps(j399))

    def test_mass_rule_matches_the_archived_components(self):
        endpoints, stats = deficit_endpoints()
        self.assertEqual(stats, {"component_count": 1230, "bands": 9833})
        self.assertEqual(len(endpoints), 11060)
        self.assertEqual(sum(weight for _, weight in endpoints), 0)

    def test_diagonal_and_wronskian_against_arb_integral(self):
        lo, hi = Fraction(1, 2), Fraction(5, 4)
        with ctx.workprec(256):
            width = arb(WIDTH.numerator) / WIDTH.denominator
            result = band_gram([(lo, -1), (hi, 1)], 2, batch_nodes=1)
            pi = arb.pi()

            def j_values(t):
                x = width * t
                first = x.sin() / x
                second = (first - x.cos()) / x
                return first, second, 3 * second / x - first

            for n, m in [(0, 0), (1, 1), (0, 2), (2, 2)]:
                def integrand(t, analytic):
                    j = j_values(t)
                    sign = (-1)**(abs(n-m)//2)
                    return 2 * width / pi * arb((2*n+1)*(2*m+1)).sqrt() * sign * j[n] * j[m]

                direct = acb.integral(integrand, acb(1) / 2, acb(5) / 4, rel_tol=arb("1e-65"))
                self.assertTrue(direct.imag.contains(0))
                self.assertTrue(result[n % 2][n // 2, m // 2].overlaps(direct.real))


if __name__ == "__main__":
    unittest.main()
