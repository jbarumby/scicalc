#!/usr/bin/env python3
"""Small, deterministic tests for the exact-rational ES8-P101 aggregator."""
import json,tempfile,unittest
from pathlib import Path
from assemble_2x2_exact import assemble

def record(kind,lo,hi):
    return {
      "status":"FINITE_INTERVAL_ENTRY_CONDITIONAL_ON_FORMULA",
      "entry_name":kind,
      "precision_bits":128,
      "epsilon_denominator":2048,
      "tolerance":"1e-8",
      "normalization":"P92",
      "entry":{"W":{"lo":lo,"hi":hi,"im_lo":"0","im_hi":"0"}}
    }
class Tests(unittest.TestCase):
    def solve(self, a, b):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/"d.json").write_text(json.dumps(a))
            (p/"c.json").write_text(json.dumps(b))
            return assemble([p/"d.json",p/"c.json"])
    def test_positive_case(self):
        d=self.solve(record("W11","44/100000","45/100000"),
                     record("W12","-30/100000","-29/100000"))
        self.assertTrue(d["determinant_positive"])
        self.assertEqual(d["status"],"INTERVAL_POSITIVE_DEFINITE_CONDITIONAL")
    def test_p99_like_indecision(self):
        d=self.solve(record("W11","256/1000000","549/1000000"),
                     record("W12","-497/1000000","-213/1000000"))
        self.assertFalse(d["determinant_positive"])
    def test_nonpositive_diagonal(self):
        d=self.solve(record("W11","-1/1000","1/1000"),
                     record("W12","-1/10000","0"))
        self.assertFalse(d["determinant_positive"])
        self.assertIsNone(d["determinant_lower"])
    def test_missing(self):
        with self.assertRaises(ValueError):
            self.solve(record("W11","1","2"),record("W11","1","2"))
    def test_bad_status(self):
        b=record("W12","-1/4","0")
        b["status"]="INCONCLUSIVE_RESOURCE_TIMEOUT"
        with self.assertRaises(RuntimeError):
            self.solve(record("W11","1","2"),b)
    def test_wrong_normalization(self):
        b=record("W12","-1/4","0")
        b["normalization"]="P89"
        with self.assertRaises(ValueError):
            self.solve(record("W11","1","2"),b)
    def test_inverted(self):
        with self.assertRaises(AssertionError):
            self.solve(record("W11","2","1"),record("W12","0","1"))
    def test_nonreal(self):
        b=record("W12","-1/4","0")
        b["entry"]["W"]["im_lo"]="1/100"
        with self.assertRaises(AssertionError):
            self.solve(record("W11","1","2"),b)
if __name__=="__main__":
    unittest.main(verbosity=2)
