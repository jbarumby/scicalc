"""Tests with real Arb enclosures; never rely on disabled `assert` under -O."""
import unittest
from fractions import Fraction

from .table import TABLE, validate_table
from .audit_b23 import audit, calc_b23, constants, q, H
from flint import ctx


class TestP76(unittest.TestCase):
    def test_table_exact_transcription_structure(self):
        self.assertTrue(validate_table())
        self.assertEqual(TABLE[0].sigma_start,16)
        self.assertEqual(TABLE[-1].sigma_start,31)
        self.assertEqual(TABLE[12].d,"0.29")
        self.assertEqual(TABLE[15].b3,"86.84")
        self.assertEqual(TABLE[0].b1,"5.069")

    def test_exact_decimal_parsing(self):
        with ctx.workprec(192):
            x=q("0.28")
            self.assertEqual(Fraction("0.28"), Fraction(7,25))
            self.assertTrue(q("0.28") > q("0.279"))
            self.assertTrue(q("0.28") < q("0.281"))

    def test_theorem_necessary_conditions(self):
        with ctx.workprec(192):
            Hball=q(str(H))
            self.assertTrue(Hball >= 1002)
            self.assertTrue(Hball >= 3*10**10)
            for row in TABLE:
                self.assertTrue(q(row.d)>0)
                self.assertTrue(q(row.delta)>0)
                self.assertTrue(q(row.kappa)>0)
                self.assertTrue(q(row.delta) < Hball.log()/20)

    def test_c5_c6_branches(self):
        with ctx.workprec(192):
            values=constants()
            self.assertTrue(values["c5"]>0)
            self.assertTrue(values["c6"]>200)
            self.assertTrue(values["m0"] < q("1.31541"))

    def test_all_16_published_b2_b3_outward(self):
        record=audit(192)
        self.assertEqual(len(record["rows"]),16)
        self.assertTrue(record["b2_b3_strict_upper_verified"])
        self.assertTrue(all(row["both_strictly_upper"] for row in record["rows"]))
        self.assertFalse(record["b1_recomputed"])
        self.assertFalse(record["full_corollary_independently_audited"])

    def test_raise_when_published_coefficient_too_small(self):
        from dataclasses import replace
        rows=list(TABLE)
        rows[0]=replace(rows[0],b2="9")
        with self.assertRaises(ArithmeticError): audit(192,tuple(rows))
        rows=list(TABLE)
        rows[15]=replace(rows[15],b3="86")
        with self.assertRaises(ArithmeticError): audit(192,tuple(rows))

    def test_256bits_strict(self):
        rec=audit(256)
        self.assertTrue(rec["b2_b3_strict_upper_verified"])
        self.assertEqual(rec["precision_bits"],256)

if __name__=="__main__":unittest.main()
