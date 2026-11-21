"""Fail-closed exact arithmetic tests for the P259 high-frequency witness gate."""
from fractions import Fraction as Q
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from .p259_witness import (
    pi_interval, bessel_polys, sqrt_lower, floor_gate, read_p258,
    prove_highpass, audit, PINNED, L
)

# The P258 artifact is downloaded by the workflow and authenticated by SHA-256.
SOURCE = Path("replay/P258_GL60.json")

class P259Tests(TestCase):
    def test_machin_enclosure(self):
        lo,hi=pi_interval()
        self.assertTrue(Q(314,100)<lo<hi<Q(315,100))
        self.assertLess(hi-lo,Q(1,10**18))

    def test_bessel_integer_polynomials(self):
        p,q=bessel_polys(3)
        self.assertEqual(p[0],[1])
        self.assertEqual(p[1],[0,1])
        self.assertEqual(q[1],[-1])
        self.assertEqual(p[2],[-1,0,3])
        self.assertEqual(q[2],[0,-3])
        self.assertEqual(p[3],[0,-6,0,15])
        self.assertEqual(q[3],[1,0,-15])

    def test_exact_sqrt_enclosures(self):
        for n in (3,7,57,115):
            lo=sqrt_lower(n,100)
            self.assertTrue(lo*lo<=n<(lo+Q(1,2**100))**2)

    def test_floor_850_strictly_valid(self):
        self.assertEqual(floor_gate(),str(Q(346,677)))
        self.assertGreater(Q(floor_gate()),Q(1,2))

    def test_pinned_sha(self):
        self.assertEqual(len(PINNED),64)
        self.assertEqual(L,Q("0.972955074527656652552676371721590"))

    def test_missing_input_fails(self):
        with TemporaryDirectory() as folder:
            with self.assertRaises(FileNotFoundError):
                read_p258(Path(folder)/"missing.json")

    def test_artifact_replay(self):
        if not SOURCE.exists():
            self.skipTest("P258 artifact absent; download Actions 44858303395")
        witness,lo,hi=read_p258(SOURCE)
        self.assertEqual(len(witness),29)
        self.assertTrue(-Q(2,10**24)<lo<hi<0)
        mass,details,upper=prove_highpass(witness)
        self.assertEqual(len(details),23)
        self.assertLess(upper,1000)
        self.assertGreater(mass,Q(95,10**25))
        self.assertGreater(lo+mass/4,Q(6,10**25))
        report=audit(SOURCE)
        self.assertFalse(report["full_377_by_377_sector_certified_positive"])
        self.assertFalse(report["Riemann_hypothesis_proved"])

    def test_tampered_artifact_rejected(self):
        if not SOURCE.exists():
            self.skipTest("P258 artifact absent")
        with TemporaryDirectory() as folder:
            p=Path(folder)/"tampered.json"
            p.write_bytes(SOURCE.read_bytes()+b"tamper")
            with self.assertRaisesRegex(ArithmeticError,"SHA-256"):
                read_p258(p)

if __name__=="__main__":
    main()
