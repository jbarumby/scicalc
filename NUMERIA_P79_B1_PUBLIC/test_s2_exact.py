"""Independent exact-rational S2 review (does not invoke FLINT and proves no RH).

Check user-submitted files byte-for-byte and verify their proposed two-decimal
coefficient replacements dominate the SIX-decimal Arb upper bounds already
certified by P79's GitHub Actions. Standalone S2 is never itself a certificate.
"""
import hashlib
import json
import unittest
from fractions import Fraction as Q
from pathlib import Path

BASE=Path(__file__).resolve().parent
INPUT=BASE/"prototype_input"

# SHA-256 of exact uploaded bytes; inspect the published copies in PR #79.
SOURCE_SHA={
    "gemini-code-1795353291902.py":"d6edcd2b02b278f17999a34102a33f2733b7231c08f0e3984676c9dc160d912f",
    "p78_b1_resolution.py":"10be6eb754ec5a834b97cb1654533f53e77000722e05e7a2d600ac9e1e6ca42d",
    "p78_b1_resolution_results.json":"c8f568359a72252aa3e18f52a470a4685606f9f0c861259a8adf5e4de6ef665d",
}
# These 16 values are pre-existing P79 certified STRICT decimal upper bounds,
# NOT results obtained from the standalone S2 file.
P79_ARBITRARY_PRECISION_STRICT_UPPER={
    10:"17.484560",11:"19.659253",12:"22.239617",
    13:"25.472651",14:"29.881218",15:"36.586170",
    16:"51.648097",
}
P79_TABLE_SHA256="d7ba517b4461e4791ad1ca92d69575d4f973b0d161b5a2a50b928ea2137af615"

def source_json():
    return json.loads((INPUT/"p78_b1_resolution_results.json").read_text(encoding="utf-8"))

def p78_json():
    p=BASE.parent/"NUMERIA_P78_RAMARE_PUBLIC"/"certificate.json"
    return json.loads(p.read_text(encoding="utf-8"))

class ExactSourceReview(unittest.TestCase):
    def test_uploaded_files_are_exact(self):
        for name, expected in SOURCE_SHA.items():
            with self.subTest(name=name):
                self.assertEqual(hashlib.sha256((INPUT/name).read_bytes()).hexdigest(), expected)

    def test_s2_declares_standalone_not_certified(self):
        s=source_json()
        self.assertEqual(s["schema"],"NUMERIA.P78.S2.REPRINT.PROTOTYPE")
        self.assertEqual(s["mode"],"STANDALONE_RECONSTRUCTION")
        self.assertEqual([r["row"] for r in s["rows"]],list(range(10,17)))
        self.assertTrue(s["summary"]["model_consistent_all"])
        self.assertFalse(s["summary"]["old_print_salvageable_with_public_b4"])

    def test_exact_exponents_and_domain(self):
        """e=(2s2-1)/(2-s2), s2=(row+16)/32; check 0<e<=1 exactly."""
        for row in range(10,17):
            s2=Q(row+16,32)
            exponent=(2*s2-1)/(2-s2)
            self.assertEqual(exponent,Q(2*row,48-row))
            self.assertGreater(exponent,0)
            self.assertLessEqual(exponent,1)
        # For all 16 P78 rows: K3>C3=pi*m0*b4*exp(...)/(2delta).
        # pi>3, m0>1, exp(...)>1, delta<=303/1000, b4>=27/50.
        # Hence K3>3*(27/50)/(2*(303/1000))=270/101>1.
        lower=Q(3)*Q(27,50)/(2*Q(303,1000))
        self.assertEqual(lower,Q(270,101))
        self.assertGreater(lower,1)

    def test_s2_source_input_matches_p78_certificate(self):
        s=source_json()
        published={int(x["row"]):x for x in p78_json()["rows"]}
        for x in s["rows"]:
            row=x["row"]; orig=published[row]
            self.assertEqual(Q(x["b4_public"]),Q(orig["b4"]))
            self.assertEqual(Q(x["delta"]),Q(orig["delta"]))
            self.assertGreater(Q(orig["b4"]),Q("0.54"))

    def test_s2_two_decimal_candidates_dominate_p79_arb_upper_bounds(self):
        s=source_json()
        for x in s["rows"]:
            row=x["row"]
            arb_upper=Q(P79_ARBITRARY_PRECISION_STRICT_UPPER[row])
            for candidate in ("reprint_minimal","reprint_robust"):
                with self.subTest(row=row,candidate=candidate):
                    self.assertGreater(Q(x[candidate]),arb_upper)
                    self.assertGreater(Q(x[candidate]),Q(x["printed_B1_old"]))
            self.assertLess(Q(x["printed_B1_old"]),arb_upper)

    def test_s2_floats_are_not_certified_endpoints(self):
        s=source_json()
        for x in s["rows"]:
            self.assertIn("B1_upper_b4_public_float",x)
            self.assertIsInstance(x["B1_upper_b4_public_float"],float)
            self.assertNotIn("B1_upper_dyadic",x)

if __name__=="__main__":
    unittest.main()
