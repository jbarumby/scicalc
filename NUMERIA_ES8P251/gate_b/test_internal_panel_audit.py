"""Reproducibility and adversarial tests for the P251 internal-panel pilot."""

from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from .internal_panel_audit import (
    PANEL,
    P245_ROOTS_SHA,
    authenticated_root_bytes,
    build_report,
    verify_report,
)

HERE = Path(__file__).resolve().parent
ARTIFACT = HERE / "data" / "P251_INTERNAL_PANEL_200_201.json"


class InternalPanelAuditTests(unittest.TestCase):
    def test_archived_panel_replay_is_byte_deterministic(self):
        archived_bytes = ARTIFACT.read_bytes()
        regenerated = (json.dumps(build_report(), indent=2, ensure_ascii=False) + "\n").encode()
        self.assertEqual(regenerated, archived_bytes)

    def test_certificate_has_six_overlapping_entries_and_explicit_scope(self):
        report = json.loads(ARTIFACT.read_text())
        verify_report(report)
        self.assertEqual(report["panel_range"], [PANEL, PANEL + 1])
        self.assertEqual(report["root_brackets_sha256"], P245_ROOTS_SHA)
        self.assertEqual(sum(len(sector["entries"]) for sector in report["sectors"]), 6)
        self.assertTrue(report["all_selected_intervals_overlap"])
        self.assertFalse(report["full_259_mode_head_recomputed"])
        self.assertFalse(report["finite_head_positivity_recertified"])
        self.assertFalse(report["Riemann_hypothesis_proved"])

    def test_wrong_panel_and_adulterated_root_brackets_fail_closed(self):
        with self.assertRaises(ValueError):
            build_report(panel=PANEL - 1)
        original = authenticated_root_bytes()
        with TemporaryDirectory() as folder:
            tampered = Path(folder) / "roots.json"
            tampered.write_bytes(original + b" ")
            with self.assertRaises(ArithmeticError):
                authenticated_root_bytes(tampered)

    def test_hash_index_and_nonoverlap_tampering_fail_closed(self):
        original = json.loads(ARTIFACT.read_text())
        bad_hash = deepcopy(original)
        bad_hash["root_brackets_sha256"] = "0" * 64
        with self.assertRaises(ArithmeticError):
            verify_report(bad_hash)

        bad_panel = deepcopy(original)
        bad_panel["panel_range"] = [199, 200]
        with self.assertRaises(ValueError):
            verify_report(bad_panel)

        bad_index = deepcopy(original)
        bad_index["sectors"][0]["entries"][0]["pair"] = [1, 1]
        with self.assertRaises(ValueError):
            verify_report(bad_index)

        bad_interval = deepcopy(original)
        bad_interval["sectors"][0]["entries"][0]["production_dyadic_bounds"] = [["1", "0"], ["1", "0"]]
        with self.assertRaises(ArithmeticError):
            verify_report(bad_interval)


if __name__ == "__main__":
    unittest.main()
