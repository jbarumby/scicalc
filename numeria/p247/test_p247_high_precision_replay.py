#!/usr/bin/env python3
"""Independent structural/rational checks for archived P247 512-bit replays."""
from __future__ import annotations

import hashlib
import json
import unittest
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
BLOCKS = ((25, 50), (500, 525), (975, 1000))


def read_report(start: int, stop: int, bits: int):
    return json.loads((HERE / "data" / f"P247D_REPLAY_{start}_{stop}_{bits}.json").read_text())


def dyadic(m: str, e: str) -> F:
    m, e = int(m), int(e)
    return F(m * (2**e)) if e >= 0 else F(m, 2**(-e))


def source_interval(ball) -> tuple[F, F]:
    (m, e), (r, k) = ball
    center, radius = dyadic(m, e), dyadic(r, k)
    return center - radius, center + radius


def interval_width(record) -> F:
    enc = record["replay_interval"]
    return F(int(enc["upper"]) - int(enc["lower"]), 2 ** enc["scale_power_of_two"])


class P247HighPrecisionReplayTests(unittest.TestCase):
    def test_all_three_artifacts_have_correct_scope_and_precision(self):
        for start, stop in BLOCKS:
            with self.subTest(block=(start, stop)):
                d = read_report(start, stop, 512)
                self.assertEqual(d["schema"], "NUMERIA.ES8P245D.SELECTED_ADDITIONAL_P237_CHUNK_REPLAY.1")
                self.assertEqual(d["parameters"]["panel_range"], [start, stop])
                self.assertEqual(d["parameters"]["fixed_point_bits"], 512)
                self.assertEqual(d["parameters"]["panels"], 25)
                self.assertEqual(d["verified_chunk_sum_enclosures"], 26)
                self.assertFalse(d["full_377_by_377_triangle_recomputed"])
                self.assertFalse(d["Weil_Zhu_identity_verified"])
                self.assertFalse(d["RH_proved"])

    def test_input_sha256_and_source_ball_containment_for_every_block(self):
        for start, stop in BLOCKS:
            d = read_report(start, stop, 512)
            for sector, source in d["reference_chunks"].items():
                with self.subTest(block=(start, stop), sector=sector):
                    path = HERE / "references" / source["filename"]
                    self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["sha256"])
                    matching = [c for c in d["comparisons"] if c["sector"] == sector]
                    self.assertEqual(len(matching), 1)
                    for comparison in matching:
                        self.assertEqual(comparison["source"]["sha256"], source["sha256"])
                        self.assertEqual(len(comparison["selected_entry_comparisons"]), 13)
                        for record in comparison["selected_entry_comparisons"]:
                            lo, hi = source_interval(record["source_Arb_ball_dyadic"])
                            self.assertEqual(str(lo), record["source_interval_lower_fraction"])
                            self.assertEqual(str(hi), record["source_interval_upper_fraction"])
                            self.assertEqual(record["replay_interval"]["scale_power_of_two"], 512)
                            self.assertLessEqual(F(int(record["replay_interval"]["lower"]), 2**512), lo)
                            self.assertGreaterEqual(F(int(record["replay_interval"]["upper"]), 2**512), hi)

    def test_high_mode_boxes_are_narrower_in_each_block(self):
        for start, stop in BLOCKS:
            old = json.loads((HERE.parent / "p245" / "gate_c" / f"P245D_SELECTED_{start}_{stop}.json").read_text())
            new = read_report(start, stop, 512)
            old_records = [r for c in old["comparisons"] for r in c["selected_entry_comparisons"] if 376 in r["pair"]]
            new_records = [r for c in new["comparisons"] for r in c["selected_entry_comparisons"] if 376 in r["pair"]]
            self.assertEqual([r["pair"] for r in old_records], [r["pair"] for r in new_records])
            self.assertEqual(len(new_records), 6)
            for old_record, new_record in zip(old_records, new_records):
                ratio = interval_width(old_record) / interval_width(new_record)
                with self.subTest(block=(start, stop), pair=new_record["pair"]):
                    self.assertGreater(ratio, 2**100)
                    self.assertLess(ratio, 2**257)


if __name__ == "__main__":
    unittest.main()
