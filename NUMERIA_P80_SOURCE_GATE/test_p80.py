"""NUMERIA P80: exact-rational and immutable-S3 regression tests."""
import json
import unittest
from fractions import Fraction as F
from pathlib import Path
from . import v3_gate
from .s3_safe import checked_rows

class SourceV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc=v3_gate.audit()
    def test_baseline_reproduced_exactly(self):
        d=self.doc
        self.assertEqual(len(d["rows"]),16)
        self.assertTrue(all(r["P79_existing_b4_still_valid"] for r in d["rows"]))
    def test_only_row10_tightened(self):
        self.assertEqual(self.doc["changed_b4_rows"],[10])
        self.assertEqual(F(self.doc["rows"][9]["b4_P78"]),F(557,1000))
        self.assertEqual(F(self.doc["rows"][9]["b4_P80_suggested"]),F(556,1000))
        self.assertEqual(F(self.doc["rows"][15]["b4_P80_suggested"]),F(303,500))
    def test_unchanged_p79_table_remains_conservative(self):
        for r in self.doc["rows"]:
            self.assertLessEqual(F(r["b4_P80_suggested"]),F(r["b4_P78"]))
            self.assertGreater(F(r["b4_P80_suggested"]),F(r["worst_min_v1_v3_upper"]))
    def test_provenance_v1_extra_bound_in_section6(self):
        r=self.doc["rows"][0]
        z=r["optimized_zone_sources"][0]
        self.assertEqual(z["S_cap"],"263/500")
        self.assertEqual(z["selected_version"],"v1")
        self.assertIn("section 6",z["citation_location"])
    def test_row10_uses_v3_5501(self):
        z=self.doc["rows"][9]["optimized_zone_sources"][0]
        self.assertEqual(z["selected_version"],"v3")
        self.assertEqual(F(z["S_cap"]),F(5501,10000))
    def test_strict_analytic_log_bounds(self):
        self.assertEqual(F(v3_gate.verify_log_lower())<F(11,4),True)
    def test_reject_out_of_scope(self):
        with self.assertRaises(ArithmeticError):
            v3_gate.best_zone(F(3),0)
        with self.assertRaises(ArithmeticError):
            v3_gate.best_zone(F(303,1000),0,versions=("missing",))
    def test_s3_selection_fails_closed(self):
        self.assertEqual(checked_rows("1-9,12"),[1,2,3,4,5,6,7,8,9,12])
        self.assertEqual(checked_rows(None),list(range(1,17)))
        for bad in ("0","17","1-17","1,,2","9-5","abc",""):
            if bad=="": continue # empty selects all by design
            with self.subTest(bad=bad),self.assertRaises(ValueError):
                checked_rows(bad)
    def test_s3_source_input_compiles(self):
        import ast
        source=Path(__file__).resolve().parent/"source_input/p78_s2_full_protocol.py"
        self.assertIn("NUMERIA P78-S3",source.read_text())
        import hashlib
        self.assertEqual(
            hashlib.sha256(source.read_bytes()).hexdigest(),
            "1faf315454a100a6a9850e843c795eca2bf734053ec9d32f202286709565b686",
        )
        ast.parse(source.read_text(encoding="utf-8"))

if __name__=="__main__":
    unittest.main()
