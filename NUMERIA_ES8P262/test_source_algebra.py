"""Exact algebra tests. No floating point, no network, no flint dependency."""
import unittest
from fractions import Fraction as Q
from pathlib import Path
from tempfile import TemporaryDirectory
import json
from .source_algebra import (
    decode_dyadic, interval, add, subtract, product, quotient, scale,
    eq417_rhs, verify
)

class SourceAlgebraTests(unittest.TestCase):
    def test_dyadic_exact(self):
        self.assertEqual(decode_dyadic(["3",-2]),Q(3,4))
        self.assertEqual(decode_dyadic(["-7",1]),-14)
    def test_intervals(self):
        self.assertEqual(add(interval(1,2),interval(3,4)),(Q(4),Q(6)))
        self.assertEqual(subtract(interval(1,2),interval(3,4)),(Q(-3),Q(-1)))
        self.assertEqual(product(interval(-1,2),interval(3,4)),(Q(-4),Q(8)))
        self.assertEqual(quotient(interval(1,2),interval(2,4)),(Q(1,4),Q(1)))
    def test_reject_denominator_through_zero(self):
        with self.assertRaises(ArithmeticError):
            quotient(interval(1,2),interval(-1,1))
    def test_eq417_symbolic_coefficient(self):
        F={0:interval(1,1),2:interval(1,1),4:interval(1,1)}
        self.assertEqual(eq417_rhs(1,F),(Q(31,4),Q(31,4)))
        self.assertEqual(eq417_rhs(2,F),(Q(191,96),Q(191,96)))
    def test_ratio_versus_printed_sign(self):
        F={0:interval(1,1),2:interval(1,1),4:interval(1,1)}
        a1=eq417_rhs(1,F);a2=eq417_rhs(2,F)
        self.assertGreater(quotient(a1,a2)[0],0)
        self.assertLess(quotient(scale(interval(63,63),Q(-1,16)),interval(383,383))[1],0)
    def test_eq417_excludes_index_zero(self):
        with self.assertRaises(ArithmeticError):
            eq417_rhs(0,{0:interval(1,1)})
    def test_wrong_source_sha_rejected(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"source.json"
            p.write_text(json.dumps({"schema":"wrong"}))
            with self.assertRaisesRegex(ArithmeticError,"SHA"):
                verify(p,expected_sha="0"*64)
    def test_wrong_source_schema_rejected(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"source.json"
            p.write_text(json.dumps({"schema":"wrong"}))
            with self.assertRaisesRegex(ArithmeticError,"schema"):
                verify(p)
    def test_bad_interval_order_rejected(self):
        with self.assertRaises(ArithmeticError):interval(2,1)
    def test_no_source_means_no_claim(self):
        with self.assertRaises(FileNotFoundError):
            verify(Path("/missing/P262_SOME_MOMENTS.json"))

if __name__=="__main__":unittest.main()
