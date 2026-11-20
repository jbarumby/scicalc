import shutil
import tempfile
import unittest
from pathlib import Path

from .audit_dependency_chain import SOURCE_HASHES, audit

ROOT = Path(__file__).parent
P255 = ROOT.parent / "NUMERIA_ES8P255" / "fixtures"


class DependencyChainTests(unittest.TestCase):
    def test_all_seven_source_files_are_pinned(self):
        result = audit(ROOT / "lean_sources", P255)
        self.assertEqual(set(result["source_sha256_verified"]), set(SOURCE_HASHES))

    def test_reduction_chain_is_present(self):
        result = audit(ROOT / "lean_sources", P255)
        self.assertIn("allChi_entry", result["checked_links"])
        self.assertIn("RH_reduction", result["checked_links"])
        self.assertIn("galoisLocalGlobalStrong_recParameter", result["concentrated_upstream_input"])

    def test_does_not_claim_kernel_replay_or_proof(self):
        result = audit(ROOT / "lean_sources", P255)
        self.assertFalse(result["lean_kernel_invoked"])
        self.assertFalse(result["full_solution_replayed"])
        self.assertFalse(result["imported_theorem_bodies_independently_verified"])

    def test_modified_source_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            sources = tmp / "src"
            shutil.copytree(ROOT / "lean_sources", sources)
            p255 = tmp / "p255"
            shutil.copytree(P255, p255)
            path = sources / "MainConclusion.lean"
            path.write_bytes(path.read_bytes() + b"\n-- altered\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                audit(sources, p255)

    def test_changed_p255_entry_point_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            p255 = tmp / "p255"
            shutil.copytree(P255, p255)
            path = p255 / "MainTheorems.lean"
            raw = path.read_bytes().replace(b"rh_of_every allChi", b"rh_of_every fakeChi")
            path.write_bytes(raw)
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                audit(ROOT / "lean_sources", p255)

    def test_missing_source_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            sources = tmp / "src"
            shutil.copytree(ROOT / "lean_sources", sources)
            (sources / "VwinOfRow.lean").unlink()
            with self.assertRaises(FileNotFoundError):
                audit(sources, P255)


if __name__ == "__main__":
    unittest.main()
