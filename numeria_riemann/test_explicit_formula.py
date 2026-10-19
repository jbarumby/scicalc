import math
import unittest

from numeria_riemann.explicit_formula import (
    first_zero_ordinates,
    psi_exact,
    psi_from_zeros_truncated,
    von_mangoldt,
)


class ExplicitFormulaTests(unittest.TestCase):
    def test_von_mangoldt_prime_powers(self):
        self.assertAlmostEqual(von_mangoldt(2), math.log(2))
        self.assertAlmostEqual(von_mangoldt(8), math.log(2))
        self.assertAlmostEqual(von_mangoldt(9), math.log(3))
        self.assertEqual(von_mangoldt(12), 0.0)
        self.assertEqual(von_mangoldt(1), 0.0)

    def test_exact_psi_inclusive(self):
        self.assertAlmostEqual(psi_exact(2), math.log(2))
        self.assertAlmostEqual(psi_exact(8), 3 * math.log(2) + math.log(3) + math.log(5) + math.log(7))
        self.assertAlmostEqual(
            psi_exact(10), 3 * math.log(2) + 2 * math.log(3) + math.log(5) + math.log(7)
        )
        self.assertAlmostEqual(psi_exact(1), 0.0)

    def test_first_zero_is_positive(self):
        zeros = first_zero_ordinates(2, dps=20)
        self.assertEqual(len(zeros), 2)
        self.assertTrue(zeros[0] > 14)
        self.assertTrue(zeros[1] > zeros[0])

    def test_truncated_formula_returns_finite_diagnostic(self):
        zeros = first_zero_ordinates(3, dps=20)
        value = psi_from_zeros_truncated(10.5, zeros, dps=20)
        self.assertTrue(math.isfinite(value))

    def test_explicit_formula_complex_zeros_do_not_double_count_pairs(self):
        x = 10.5
        ordinate_form = psi_from_zeros_truncated(
            x, [14.134725141734695], dps=30, assume_critical_line=True
        )
        complex_pair_form = psi_from_zeros_truncated(
            x,
            [complex(0.5, 14.134725141734695), complex(0.5, -14.134725141734695)],
            dps=30,
            assume_critical_line=False,
        )
        self.assertAlmostEqual(ordinate_form, complex_pair_form, places=12)


if __name__ == "__main__":
    unittest.main()
