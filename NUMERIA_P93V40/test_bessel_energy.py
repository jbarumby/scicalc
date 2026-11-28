"""Independent Arb comparisons for the energy and ratio recurrences."""

import unittest

from flint import arb, ctx

from .bessel_energy import bessels_at_point


class EnergyTests(unittest.TestCase):
    def test_high_degree_against_direct_bessel(self):
        with ctx.workprec(384):
            for x, last, check in (
                (arb(200), 399, (0, 1, 199, 200, 201, 399)),
                (arb(4781) / 4, 1200, (0, 1, 1194, 1195, 1196, 1200)),
                (arb(4097), 1200, (0, 1, 1198, 1199, 1200)),
            ):
                with self.subTest(x=x, last=last):
                    actual = bessels_at_point(x, last)
                    self.assertEqual(len(actual), last + 1)
                    self.assertTrue(all(v.rad() < arb(2)**-330 for v in actual))
                    with ctx.workprec(2048):
                        scale = (arb.pi() / (2 * x)).sqrt()
                        for n in check:
                            direct = scale * x.bessel_j(arb(2 * n + 1) / 2)
                            self.assertTrue(direct.rad() < arb(2)**-400)
                            self.assertTrue(actual[n].overlaps(direct))

    def test_small_argument_and_endpoint(self):
        with ctx.workprec(256):
            x = arb(1) / 4
            values = bessels_at_point(x, 8)
            self.assertTrue(bessels_at_point(x, 0)[0].overlaps(x.sin() / x))
            self.assertTrue(all(v > 0 for v in values))
            for a, b in zip(values, values[1:]):
                self.assertTrue(a > b)

    def test_invalid_inputs(self):
        for x, degree in ((arb(0), 0), (arb(1, 1), 0), (arb(1), -1), (arb(1), True)):
            with self.subTest(x=x, degree=degree):
                with self.assertRaises(ValueError):
                    bessels_at_point(x, degree)


if __name__ == "__main__":
    unittest.main()
