#!/usr/bin/env python3
"""Small deterministic checks for the P245D chunk replay harness."""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from p245d_additional_chunk_replay import (
    N, UPPER, PINNED_HASHES, SUPPORTED_STARTS, high_mode_bound, triangle_index,
)
from p245b_certified_gl60 import I, PREC, S, source_ball_to_fraction


class P245DReplayTests(unittest.TestCase):
    def test_upper_triangle_index_bijection_endpoints(self):
        self.assertEqual(triangle_index(0, 0), 0)
        self.assertEqual(triangle_index(0, 1), 1)
        self.assertEqual(triangle_index(1, 1), N)
        self.assertEqual(triangle_index(376, 376), UPPER - 1)

    def test_exactly_three_supported_blocks_and_six_pins(self):
        self.assertEqual(SUPPORTED_STARTS, (25, 500, 975))
        self.assertEqual(len(PINNED_HASHES), 6)
        self.assertTrue(all(len(value) == 64 for value in PINNED_HASHES.values()))

    def test_high_mode_bound_is_positive_and_ordered(self):
        for parity in (0, 1):
            bounds = [high_mode_bound(panel, parity) for panel in (25, 500, 975)]
            self.assertTrue(all(bound > 0 for bound in bounds))
            self.assertLess(bounds[0], bounds[1])
            self.assertLess(bounds[1], bounds[2])

    def test_archived_reference_balls_are_enclosed(self):
        for start in SUPPORTED_STARTS:
            stop = start + 25
            path = HERE / f"P245D_SELECTED_{start}_{stop}.json"
            data = json.loads(path.read_text())
            self.assertEqual(data["parameters"]["panel_range"], [start, stop])
            self.assertEqual(data["verified_chunk_sum_enclosures"], 26)
            for comparison in data["comparisons"]:
                sector = comparison["sector"]
                self.assertEqual(comparison["source"]["sha256"], PINNED_HASHES[start, sector])
                self.assertEqual(len(comparison["selected_entry_comparisons"]), 13)
                for record in comparison["selected_entry_comparisons"]:
                    lo, hi = source_ball_to_fraction(record["source_Arb_ball_dyadic"])
                    self.assertEqual(str(lo), record["source_interval_lower_fraction"])
                    self.assertEqual(str(hi), record["source_interval_upper_fraction"])
                    enclosure = record["replay_interval"]
                    self.assertEqual(enclosure["scale_power_of_two"], PREC)
                    replay = I(int(enclosure["lower"]), int(enclosure["upper"]))
                    self.assertTrue(replay.contain_fraction(lo))
                    self.assertTrue(replay.contain_fraction(hi))


if __name__ == "__main__":
    unittest.main()
