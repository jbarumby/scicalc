"""P251 transfer gate regression and fail-closed checks."""
import unittest
from fractions import Fraction as F
from .p251_transfer_gate import certify, need


class P251TransferGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = certify()

    def test_nontransfer_proved_by_exact_rational_ratio(self):
        r = self.result
        self.assertGreater(F(r['H00_lower_exact']), F(r['old_margin_exact']))
        self.assertGreater(F(r['H00_to_old_margin_ratio_exact']), F(2*10**8))

    def test_four_sine_cells_are_bounded(self):
        self.assertEqual(len(self.result['four_sine_cell_upper_bounds_exact']), 4)
        self.assertTrue(all(F(bound) < -F(1, 2)
                            for bound in self.result['four_sine_cell_upper_bounds_exact']))

    def test_no_false_global_claims(self):
        self.assertFalse(self.result['reduced_259_mode_matrix_certified'])
        self.assertFalse(self.result['reduced_matrix_indefinite_proved'])
        self.assertFalse(self.result['RH_proved'])

    def test_fail_closed_works_under_python_optimized_mode(self):
        with self.assertRaises(ArithmeticError):
            need(False, 'malformed certificate')

if __name__ == '__main__':
    unittest.main()
