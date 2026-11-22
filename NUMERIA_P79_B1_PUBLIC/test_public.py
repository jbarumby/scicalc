"""Unit tests for P79 certificate inputs and exact rounding (no quadrature here)."""
import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch
from . import audit_public as p79

class P79ExactTests(unittest.TestCase):
    def test_dyadic_exact(self):
        self.assertEqual(p79.exact_from_dyadic(["3",-1]),Fraction(3,2))
        self.assertEqual(p79.exact_from_dyadic(["5",2]),Fraction(20))

    def test_round_up_strictly(self):
        self.assertEqual(p79.rounded_up_dyadic(["3",-1]),"1.500001")
        self.assertEqual(p79.rounded_up_dyadic(["7",-3]),"0.875001")
        self.assertGreater(Fraction(p79.rounded_up_dyadic(["1",-20])),Fraction(1,2**20))

    def test_source_rows(self):
        rows,sha=p79.load_source()
        self.assertEqual(len(rows),16)
        self.assertEqual(len(sha),64)
        self.assertEqual(sum(Fraction(x["b4"])==Fraction("0.54") for x in rows),9)
        self.assertEqual(rows[15]["b4"],"303/500")

    def test_tampered_source_fails(self):
        rows=json.loads(p79.P78_CERT.read_text())
        rows["rows"][9]["b4"]="27/50"
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"bad.json"
            p.write_text(json.dumps(rows))
            with patch.object(p79,"P78_CERT",p):
                with self.assertRaises(ArithmeticError):
                    p79.load_source()

    def test_source_wrong_row_fails(self):
        rows=json.loads(p79.P78_CERT.read_text())
        rows["rows"][0]["row"]=42
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"bad.json"
            p.write_text(json.dumps(rows))
            with patch.object(p79,"P78_CERT",p):
                with self.assertRaises(ArithmeticError):
                    p79.load_source()

if __name__=="__main__":
    unittest.main()
