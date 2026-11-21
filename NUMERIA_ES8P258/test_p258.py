"""Regression tests for the raised-floor P258 calculations."""

from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from flint import arb, arb_mat, ctx

from .audit import update_head, _verify_tail
from .concentration import (
    aggregate_chunks, analytic_concentration, compare_concentration, concentration_chunk,
    load_authenticated, write_chunk,
)
from .independent_rayleigh import replay_rayleigh
from .strip_budget import certify_strip_constants


class RaisedFloorTest(unittest.TestCase):
    def test_independent_concentration_and_parity(self):
        options = {"half_width": Fraction(19, 20), "cutoff": 1, "modes": 3}
        sampled = concentration_chunk(0, 2, **options)
        endpoint = analytic_concentration(**options)
        result = compare_concentration(sampled, endpoint, **options)
        self.assertEqual([result[s]["endpoint_enclosures_contained"] for s in ("even", "odd")], [6, 6])
        self.assertTrue(endpoint[0][0, 0] > endpoint[1][0, 0] > 0)
        self.assertTrue(endpoint[0][0, 0] < 1)

    def test_update_does_not_double_count_the_floor(self):
        with ctx.workprec(192):
            head = arb_mat([[1, 2], [2, 3]])
            gram = arb_mat([[arb(1) / 2, arb(1) / 4], [arb(1) / 4, arb(3) / 4]])
            updated = update_head(head, gram, arb(1) / 8)
            self.assertEqual(updated[0, 0], arb(17) / 16)
            self.assertEqual(updated[0, 1], arb(63) / 32)
            self.assertEqual(updated[1, 1], arb(97) / 32)

    def test_exact_rayleigh_sign_and_errors(self):
        rows = [
            [0, 0, [["-2", "0"], ["-1", "0"]]],
            [0, 2, [["1", "0"], ["3", "0"]]],
            [2, 2, [["-1", "0"], ["0", "0"]]],
        ]
        result = replay_rayleigh(("-1/2", 0, 2), {"assembled": (rows, 1)}, operator_errors={"GL": "1/4"})
        self.assertEqual(result["rayleigh_interval_exact"], ["-185/68", "-19/68"])
        self.assertEqual(result["witness_strict_sign"], "negative")
        with self.assertRaises(ValueError):
            replay_rayleigh(("-1/2", 0, 2), {"assembled": (rows[::2], 1)})

    def test_rejects_changed_aggregate_or_coverage_gap(self):
        with TemporaryDirectory() as work:
            shard = Path(work) / "chunk.gz"
            dest = Path(work) / "aggregate.gz"
            options = {"cutoff": 1, "modes": 1}
            write_chunk(shard, 0, 1, **options)
            with self.assertRaisesRegex(ValueError, "incomplete"):
                aggregate_chunks([shard], dest, **options)
            write_chunk(shard, 0, 2, **options)
            aggregate_chunks([shard], dest, **options)
            from hashlib import sha256
            digest = sha256(dest.read_bytes()).hexdigest()
            load_authenticated(dest, digest, **options)
            dest.write_bytes(dest.read_bytes() + b"changed")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                load_authenticated(dest, digest, **options)

    def test_tail_and_strip(self):
        self.assertEqual(len(_verify_tail()), 64)
        constants = certify_strip_constants()
        self.assertEqual(constants["strip_loss_factor_exact_upper"], "23785635/16799")
        self.assertEqual(constants["strip_loss_factor_integer_upper"], 1416)


if __name__ == "__main__":
    unittest.main()
