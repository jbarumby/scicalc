"""Exact-rational regression and explicit failure checks for P78."""
import unittest
from fractions import Fraction as F
from .restricted_lemma6 import DELTAS, certificate, elementary, row_result

class TestP78(unittest.TestCase):
    def test_log_bounds(self):
        self.assertTrue(elementary())
    def test_all_rows(self):
        data=certificate()
        self.assertEqual(len(data["rows"]),16)
        self.assertTrue(all(F(r["margin"])>0 for r in data["rows"]))
        self.assertFalse(data["B1_recomputed"])
        self.assertFalse(data["global_RH_proved"])
    def test_nine_original_bounds(self):
        rows=certificate()["rows"]
        self.assertEqual([F(r["b4"])==F(".54") for r in rows],
                         [True]*9+[False]*7)
    def test_last_row(self):
        row=certificate()["rows"][-1]
        self.assertEqual(F(row["b4"]),F(".606"))
        self.assertEqual(F(row["margin"]),F(2321,14000000))
    def test_reject_parameters_outside_scope(self):
        for delta in (F(0),F(-1),F(".304")):
            with self.assertRaises(ValueError):
                row_result(1,delta)
    def test_larger_than_private_lemma_on_last_row(self):
        row=certificate()["rows"][-1]
        self.assertLess(F(".54"),F(row["worst_upper"]))
        self.assertLess(F(".605"),F(row["worst_upper"]))
    def test_exact_rational_only(self):
        for row in certificate()["rows"]:
            self.assertTrue(all("/" in row[key] for key in ("delta","b4","margin")))

if __name__=="__main__":
    unittest.main()
