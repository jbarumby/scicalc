"""Independent normalization and rational bound checks for the wider window."""

from fractions import Fraction
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from flint import arb, arb_mat, ctx

from numeria.p244.package.p244_analytic_gates import exp_lower, exp_upper
from NUMERIA_ES8P249.legendre_tail import fixed_symbol_half_width, modal_tail_bounds
from NUMERIA_ES8P251.gate_b.finite_head import aggregate_chunks
from NUMERIA_ES8P251.gate_b.finite_head_certify import assemble_finite_head, gl60_error_budget
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import chunk_report, fixed_beta


L = Fraction(19, 20)


class WiderWindowTests(unittest.TestCase):
    def test_fixed_symbol_and_modal_budgets(self):
        self.assertGreater(exp_lower(Fraction(161, 100), 8), 5)
        self.assertLess(exp_upper(2 * L, 8), 7)
        self.assertEqual(fixed_symbol_half_width(L), L)
        with self.assertRaises(ValueError):
            fixed_symbol_half_width(Fraction(1))
        for parity, sector in enumerate(("even", "odd")):
            bound = gl60_error_budget(sector, half_width=L)
            self.assertEqual(bound["argument_upper"], 294)
            self.assertLess(bound["row_sum_bound"], Fraction(1391, 10**50))
            coupling_squared, deviation = modal_tail_bounds(309, 259, parity, half_width=L)
            self.assertLess(coupling_squared, Fraction(1, 10**108))
            self.assertLess(deviation, Fraction(1, 10**117))

    def test_wider_constant_pole_agrees_with_direct_moment_integral(self):
        with ctx.workprec(256):
            matrix, _ = assemble_finite_head(arb_mat(1, 1), "even", half_width=L)
            width = arb(L.numerator) / L.denominator
            pole = (2 * width).sqrt() * (width / 2).sinh() / (width / 2)
            self.assertTrue(matrix[0, 0].overlaps(fixed_beta() + 2 * pole**2))

    def test_aggregate_rejects_mixed_widths(self):
        with TemporaryDirectory() as directory:
            combined = []
            for panel, width in ((0, L), (1, Fraction(161, 200))):
                chunk = chunk_report(panel, panel + 1, cutoff=1, modes=1, half_width=width)
                path = Path(directory) / f"{panel}.json.gz"
                path.write_bytes(gzip.compress(json.dumps(chunk).encode(), mtime=0))
                combined.append(path)
            with self.assertRaisesRegex(ValueError, "half_width"):
                aggregate_chunks(combined, cutoff=1, modes=1, half_width=L)


if __name__ == "__main__":
    unittest.main()
