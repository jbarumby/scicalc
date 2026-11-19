#!/usr/bin/env python3
"""Regression tests for the replayable P237 analytic Gate C lemmas."""
from __future__ import annotations

import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from p237_analytic_gates import make_gates, need
from p237_gl60_analytic_certificate import build, ceil_sqrt
from p237_symbol_majorant import check_certificate, exp_lower


class P237AnalyticCertificateTests(unittest.TestCase):
    def test_fail_closed_gate_is_not_a_python_assert(self):
        with self.assertRaises(ArithmeticError):
            need(False, "deliberate failing gate")

    def test_exact_bernstein_and_argument_geometry(self):
        result = make_gates()
        self.assertTrue(all(result["checks"].values()))
        self.assertEqual(result["t_real_abs_upper"], "3025/6")
        self.assertEqual(result["t_imag_abs_upper"], "1/3")
        self.assertEqual(result["digamma_real_range"], ["1/12", "5/12"])
        self.assertEqual(result["bessel_Lt_real_abs_upper"], "19481/48")
        self.assertEqual(result["bessel_Lt_imag_abs_upper"], "161/600")
        self.assertLess(F(result["bessel_Lt_modulus_squared_upper"]), 406**2)

    def test_symbol_bound_and_beta_gate(self):
        result = check_certificate()
        self.assertEqual(result["symbol_minus_beta_bound_strict_less"], "274/7")
        self.assertTrue(result["beta_rational_lo_gt_zero"])
        self.assertTrue(result["beta_rational_hi_lt_two"])
        self.assertFalse(result["RH_proved"])
        self.assertGreater(exp_lower(6, 7), 253)

    def test_exact_square_root_ceiling(self):
        self.assertEqual(ceil_sqrt(F(16, 9)), 2)
        self.assertEqual(ceil_sqrt(F(17, 9)), 2)
        self.assertEqual(ceil_sqrt(F(37, 9)), 3)

    def test_gl60_budget_both_parities(self):
        result = build()
        self.assertEqual(result["per_panel_error_coefficient_strict_less"], "50 / 3^120")
        self.assertEqual(result["Gauss_legendre_nodes"], 60)
        self.assertEqual(result["degree_of_exactness"], 119)
        self.assertEqual(set(result["sectors"]), {"even", "odd"})
        for sector in result["sectors"].values():
            self.assertTrue(sector["row_sum_below_1e-46"])
            self.assertTrue(sector["last_mode_GL60_diagonal_error_below_1e-200"])
        self.assertFalse(result["RH_proved"])


if __name__ == "__main__":
    unittest.main()
