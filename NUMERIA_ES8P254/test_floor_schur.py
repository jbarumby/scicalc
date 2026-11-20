"""Fail-closed regression suite for P254 exact-rational infrastructure."""
import unittest
from fractions import Fraction as F
from .floor_schur import (iv, plus, minus, times, dyadic, delta_from_beta,
                          floor_upgraded_intervals, gershgorin_positive_lower, block_schur_gate)


class FloorBridge(unittest.TestCase):
    def test_interval_product_crossing_zero(self):
        self.assertEqual(times((-2, 3), (4, 5)), (F(-10), F(15)))
        self.assertEqual(minus((1, 2), (-3, 4)), (F(-3), F(5)))

    def test_beta_interval_and_raised_floor(self):
        self.assertEqual(delta_from_beta((F(1,1000), F(2,1000))), (F(248,1000),F(249,1000)))

    def test_floor_not_validated_at_309(self):
        with self.assertRaises(ValueError):
            delta_from_beta((F(1,1000),F(2,1000)),cutoff=309)

    def test_negative_delta_must_fail(self):
        with self.assertRaises(ArithmeticError):
            delta_from_beta((F(1,3), F(1,2)))

    def test_exact_two_by_two_reassembly(self):
        old = [[F(1,4),F(1,5)],[F(1,5),F(1,2)]]
        gram = [[F(3,4),F(1,10)],[F(1,10),F(1,4)]]
        new = floor_upgraded_intervals(old,gram,(F(1,4),F(1,4)))
        self.assertEqual(new[0][0],(F(5,16),F(5,16)))
        self.assertEqual(new[0][1],(F(7,40),F(7,40)))
        self.assertEqual(new[1][1],(F(11,16),F(11,16)))

    def test_all_beta_and_gram_interval_extremes(self):
        old = [[(F(1,2),F(1,2))]]
        gram = [[(F(1,4),F(1,3))]]
        new = floor_upgraded_intervals(old,gram,(F(1,10),F(1,5)))
        self.assertEqual(new[0][0],(F(17,30),F(13,20)))

    def test_asymmetric_data_must_fail(self):
        with self.assertRaises(ValueError):
            floor_upgraded_intervals([[1,2],[3,4]],[[0,0],[0,0]],(1,1))

    def test_positive_gershgorin(self):
        a = [[(2,2),(-F(1,10),F(1,10))],[(-F(1,10),F(1,10)),(3,3)]]
        self.assertEqual(gershgorin_positive_lower(a),F(19,10))

    def test_gershgorin_rejects_borderline(self):
        with self.assertRaises(ArithmeticError):
            gershgorin_positive_lower([[(1,1),(-1,1)],[(-1,1),(1,1)]])

    def test_schur_beats_naive_weyl_on_toy(self):
        # a-||B||=1/10-1/2<0, but strict block Schur succeeds.
        gate=block_schur_gate(F(1,10),F(10),[[(F(1,2),F(1,2))]],F(1,20))
        self.assertTrue(gate["certified_for_supplied_bounds"])
        self.assertGreater(F(gate["strict_schur_slack"]),0)

    def test_schur_rejects_large_cross_block(self):
        with self.assertRaises(ArithmeticError):
            block_schur_gate(F(1,10),F(10),[[(2,2)]])

    def test_schur_rejects_unproved_diagonal(self):
        with self.assertRaises(ArithmeticError):
            block_schur_gate(F(0),F(10),[[F(0)]])

    def test_dyadic_exact(self):
        self.assertEqual(dyadic(["3","-2"]),F(3,4))
        self.assertEqual(dyadic(["-1","3"]),F(-8))


if __name__ == "__main__":
    unittest.main()
