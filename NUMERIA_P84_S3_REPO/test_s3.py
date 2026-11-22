"""P84 fail-closed tests for versioned TABLE and exact source provenance."""
import unittest
from fractions import Fraction as F

from NUMERIA_P76_COEFFICIENTS.table import (
    TABLE, TABLE_CS_V2_PRINTED, TABLE_S2_PUBLIC_ROBUST,
    validate_table_s2_robust
)
from NUMERIA_P75_CS2.cs_v2 import TABLE_V2, TABLE_V2_S2_PUBLIC_ROBUST
from .audit_s3 import strict_source_checks, exact, DOCS, DOC_SHA, ORIG
from .scan_downstream import scan

ROBUST=("17.49","19.67","22.25","25.48","29.89","36.59","51.65")
PRINTED=("17.21","19.31","21.77","24.83","29.01","35.35","46.06")

class S3AdoptionTests(unittest.TestCase):
    def test_versioned_tables_and_immutability(self):
        self.assertTrue(validate_table_s2_robust())
        self.assertIs(TABLE,TABLE_CS_V2_PRINTED)
        self.assertEqual(tuple(r.b1 for r in TABLE[9:]),PRINTED)
        self.assertEqual(tuple(r.b1 for r in TABLE_S2_PUBLIC_ROBUST[9:]),ROBUST)
        self.assertEqual(tuple(row[0] for row in TABLE_V2_S2_PUBLIC_ROBUST[9:]),ROBUST)
        self.assertEqual(tuple(row[0] for row in TABLE_V2[9:]),PRINTED)

    def test_exact_source_and_patch(self):
        rows,hashval,patch=strict_source_checks()
        self.assertEqual(len(rows),16)
        self.assertEqual(len(hashval),64)
        self.assertEqual(len(patch["rows"]),7)

    def test_no_float_endpoints(self):
        self.assertEqual(exact(["3",-1]),F(3,2))
        self.assertEqual(exact(["5",2]),F(20))
        for bad in ([1,-1],["1",1.0],["1"],["1",100001]):
            with self.subTest(bad=bad),self.assertRaises((ArithmeticError,ValueError)):
                exact(bad)

    def test_legacy_reports_are_sources_not_certificates(self):
        import json
        for doc in DOCS:
            x=json.loads((ORIG/doc).read_text())
            self.assertIn("REPO_STANDARD",x["schema"])
            self.assertIn("envelopes",x["padrao"])
            self.assertTrue(isinstance(x["padrao"]["envelopes"]["alfa_beta_D_rel"],float))
            self.assertEqual(len(DOC_SHA[doc]),64)

    def test_real_repo_downstream_scan(self):
        d=scan()
        self.assertTrue(d["original_kit_scan_does_not_certify_repo"])
        self.assertEqual(d["historical_kit_scan_root"],"/tmp/kit_p78")
        self.assertGreater(d["files_examined"],17)
        self.assertGreater(d["files_with_B1_references"],0)
        self.assertIn("VERSIONED_TABLE_OR_CONSUMER",d["classifications"])

    def test_robust_exceeds_existing_p79_six_decimal_majorants(self):
        upper=("17.484560","19.659253","22.239617","25.472651",
               "29.881218","36.586170","51.648097")
        for robust,old,ub in zip(ROBUST,PRINTED,upper):
            self.assertGreater(F(robust),F(ub))
            self.assertLess(F(old),F(ub))

if __name__=="__main__":
    unittest.main()
