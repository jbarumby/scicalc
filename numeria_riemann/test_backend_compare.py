import unittest

from numeria_riemann.backend_compare import zero_residual_report


class BackendCompareTests(unittest.TestCase):
    def test_zero_residual_report_shape(self):
        rows = zero_residual_report(2, dps=20)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["index"], 1)
        self.assertGreater(rows[0]["ordinate"], 14)
        self.assertIn("not_zero_certificate", rows[0]["status"])

    def test_optional_flint_error_is_clear(self):
        try:
            import flint  # noqa: F401
        except ImportError:
            from numeria_riemann.backend_compare import compare_zeta_values
            with self.assertRaisesRegex(ImportError, "python-flint is optional"):
                compare_zeta_values([(2.0, 0.0)], dps=20)


if __name__ == "__main__":
    unittest.main()
