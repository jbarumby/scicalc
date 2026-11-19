#!/usr/bin/env python3
"""Independent structural/rational checks for the archived P247 512-bit replay."""
from __future__ import annotations

import hashlib
import json
import unittest
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARTIFACT = HERE / "data" / "P247D_REPLAY_25_50_512.json"


def dyadic(m: str, e: str) -> F:
    m, e = int(m), int(e)
    return F(m * (2**e)) if e >= 0 else F(m, 2**(-e))


def source_interval(ball) -> tuple[F, F]:
    (m, e), (r, k) = ball
    center, radius = dyadic(m, e), dyadic(r, k)
    return center - radius, center + radius


class P247HighPrecisionReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(ARTIFACT.read_text())

    def test_artifact_scope_and_precision_are_explicit(self):
        d = self.data
        self.assertEqual(d["schema"], "NUMERIA.ES8P245D.SELECTED_ADDITIONAL_P237_CHUNK_REPLAY.1")
        self.assertEqual(d["parameters"]["panel_range"], [25, 50])
        self.assertEqual(d["parameters"]["fixed_point_bits"], 512)
        self.assertEqual(d["parameters"]["panels"], 25)
        self.assertEqual(d["verified_chunk_sum_enclosures"], 26)
        self.assertFalse(d["full_377_by_377_triangle_recomputed"])
        self.assertFalse(d["Weil_Zhu_identity_verified"])
        self.assertFalse(d["RH_proved"])

    def test_input_sha256_and_source_ball_containment(self):
        d = self.data
        for sector, source in d["reference_chunks"].items():
            path = HERE / "references" / source["filename"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["sha256"])
            matching = [c for c in d["comparisons"] if c["sector"] == sector]
            self.assertTrue(matching)
            for comparison in matching:
                self.assertEqual(comparison["source"]["sha256"], source["sha256"])
                self.assertEqual(len(comparison["selected_entry_comparisons"]), 13)
                for record in comparison["selected_entry_comparisons"]:
                    lo, hi = source_interval(record["source_Arb_ball_dyadic"])
                    self.assertEqual(str(lo), record["source_interval_lower_fraction"])
                    self.assertEqual(str(hi), record["source_interval_upper_fraction"])
                    enclosure = record["replay_interval"]
                    self.assertEqual(enclosure["scale_power_of_two"], 512)
                    replay_lo = F(int(enclosure["lower"]), 2**512)
                    replay_hi = F(int(enclosure["upper"]), 2**512)
                    self.assertLessEqual(replay_lo, lo)
                    self.assertGreaterEqual(replay_hi, hi)

    def test_high_mode_rounding_boxes_shrink_by_about_2_to_256(self):
        baseline_path = HERE.parent / "p245" / "gate_c" / "P245D_SELECTED_25_50.json"
        baseline = json.loads(baseline_path.read_text())

        def widths(report):
            result = []
            for comparison in report["comparisons"]:
                for record in comparison["selected_entry_comparisons"]:
                    if 376 not in record["pair"]:
                        continue
                    enc = record["replay_interval"]
                    scale = 2 ** enc["scale_power_of_two"]
                    result.append(F(int(enc["upper"]) - int(enc["lower"]), scale))
            return result

        old, new = widths(baseline), widths(self.data)
        self.assertEqual(len(old), 6)
        self.assertEqual(len(new), 6)
        for old_width, new_width in zip(old, new):
            self.assertGreater(old_width / new_width, 2**255)
            self.assertLess(old_width / new_width, 2**257)


if __name__ == "__main__":
    unittest.main()
