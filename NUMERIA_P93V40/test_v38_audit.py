"""Independent arithmetic, tampering and transfer checks for the V.38 replay."""

from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from flint import ctx
import numpy as np

from .v38_audit import ATTACHMENT, GRID, _GateReplay, _penalties, components, run, schur_bound


class V38AuditTests(unittest.TestCase):
    def test_beta_shift_is_applied_before_grouping(self):
        rows = [(25600, 25601, 250002), (25601, 25602, 250006),
                (25602, 25603, 250000), (25603, 25604, 250003)]
        self.assertEqual(components(rows, Fraction(7, 4), 1600),
                         [(25600, 25602, 8, 6), (25603, 25604, 3, 3)])

    def test_integer_schur_replay_survives_machine_overflow(self):
        kernel = np.array([[3, 4], [4, 5]], dtype=np.int64)
        masses = [8, 12]
        weights = [2**63 + 30, 2**64 + 1]
        exact_first = Fraction(106 * (24 * weights[0] + 48 * weights[1]), 333 * GRID * 16_000_000 * weights[0])
        exact_second = Fraction(106 * (32 * weights[0] + 60 * weights[1]), 333 * GRID * 16_000_000 * weights[1])
        self.assertEqual(schur_bound(kernel, masses, weights), max(exact_first, exact_second))
        with self.assertRaisesRegex(ValueError, "positive integral"):
            schur_bound(kernel, masses, [0, 1])

    def test_a4_first_dimension_and_underflow(self):
        c, d, lam = map(Fraction, ("2791817594563/1509949440000", "3870563/1000000", "9/10"))
        _, f1, eta = _penalties(c, d, lam, 20)
        with ctx.workprec(192):
            replay = _GateReplay()
            rho = Fraction(5, 4)
            self.assertFalse(replay.check(4500, rho, 3325, 0, eta, 50)[0])
            self.assertTrue(replay.check(4500, rho, 3350, 0, eta, 50)[0])
            passes, g, b2 = replay.check(4500, rho, 4100, 1, f1, 50 + 10*d)
            self.assertTrue(passes)
            self.assertTrue(g > 0)
            self.assertTrue(b2 > 0)
            self.assertEqual(float(g), 0.0)
            self.assertEqual(float(b2), 0.0)

    def test_attachment_tampering_is_rejected_before_replay(self):
        with TemporaryDirectory() as directory:
            changed = Path(directory) / ATTACHMENT.name
            changed.write_bytes(ATTACHMENT.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                run(changed)


if __name__ == "__main__":
    unittest.main()
