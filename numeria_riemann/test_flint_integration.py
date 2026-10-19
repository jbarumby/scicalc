import unittest

from numeria_riemann.backend_compare import compare_zeta_values, compare_zero_ordinates


class FlintIntegrationTests(unittest.TestCase):
    def test_arb_zeta_and_zero_enclosures(self):
        try:
            import flint  # noqa: F401
        except ImportError:
            self.skipTest("optional python-flint dependency not installed")
        values = compare_zeta_values([(2.0, 0.0)], dps=20)
        self.assertEqual(len(values), 1)
        self.assertLess(values[0]["absolute_difference"], 1e-12)
        zeros = compare_zero_ordinates(1, dps=20)
        self.assertEqual(len(zeros), 1)
        self.assertLess(zeros[0]["midpoint_difference"], 1e-12)
        self.assertGreaterEqual(zeros[0]["flint_radius_upper_bound"], 0.0)


if __name__ == "__main__":
    unittest.main()
