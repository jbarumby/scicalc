import unittest
from fractions import Fraction as F

from flint import arb, ctx

import p246_cofinal_scaling as p246


class P246CofinalScalingTests(unittest.TestCase):
    def test_rational_upper_cutoff_strictly_dominates_arb_ball(self):
        ctx.prec = 256
        requested = 4 * arb.pi() * p246.aq(F(161, 200)).exp()
        rational_upper = p246.arb_upper_fraction(requested) + 1
        gap = p246.aq(rational_upper) - requested
        self.assertGreater(gap.lower(), 0)

    def test_p237_half_width_has_expected_prime_power_cutoff(self):
        row = p246.cofinal_row(F(161, 200), 256)
        self.assertEqual(row["prime_power_cutoff_max_n"], 5)
        self.assertEqual(row["number_of_prime_powers_in_sum"], 4)
        self.assertGreaterEqual(F(row["T_rational_upper_cutoff"]), F(504))
        self.assertLessEqual(F(161, 200) * F(row["T_rational_upper_cutoff"]),
                             F(row["modes_per_parity_sufficient"]))
        self.assertTrue(all(s["ratio_below_3_over_40"]
                            for s in row["sector_ratio_checks"]))
        self.assertTrue(all(F(s["ratio_squared_exact"]) < F(9, 1600)
                            for s in row["sector_ratio_checks"]))

    def test_larger_window_cutoff_is_exact_and_cofinal_scale_rises(self):
        small = p246.cofinal_row(F(1), 256)
        large = p246.cofinal_row(F(2), 256)
        self.assertEqual(small["prime_power_cutoff_max_n"], 7)
        self.assertEqual(large["prime_power_cutoff_max_n"], 54)
        self.assertEqual(large["number_of_prime_powers_in_sum"], 24)
        self.assertGreater(int(large["modes_per_parity_sufficient"]),
                           int(small["modes_per_parity_sufficient"]))
        for row in (small, large):
            self.assertTrue(all(F(s["ratio_squared_exact"]) < F(9, 1600)
                                for s in row["sector_ratio_checks"]))

    def test_nonpositive_width_fails_closed(self):
        with self.assertRaises(ValueError):
            p246.cofinal_row(F(0), 256)


if __name__ == "__main__":
    unittest.main()
