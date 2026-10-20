import unittest

from flint import arb

from numeria_riemann.weil_interval_audit import _omega, _sinc, run


class WeilIntervalAuditTests(unittest.TestCase):
    def test_sinc_interval_contains_zero_singularity(self):
        value = _sinc(arb(0, "0.01"))
        self.assertTrue(value.is_finite())
        self.assertTrue(value.lower() < arb("1.001"))
        self.assertTrue(value.upper() > arb("0.999"))

    def test_omega_interval_is_finite(self):
        value = _omega(arb("1 +/- 0.001"), terms=32)
        self.assertTrue(value.is_finite())

    def test_small_interval_audit_returns_explicit_status(self):
        report = run(L_value=0.45, N=4, R_value=40, step_value=0.05,
                     terms=32, precision=128)
        self.assertIn(report["status"], {
            "certified_positive_definite_finite_matrix",
            "inconclusive_interval_ldl",
        })
        self.assertEqual(report["projected_dimension"], 2)
        self.assertEqual(report["parameters"]["quadrature_cells"], 800)
        self.assertGreaterEqual(len(report["interval_ldl_pivots"]), 1)
self.assertLessEqual(len(report["interval_ldl_pivots"]), 2)
if report["status"] == "certified_positive_definite_finite_matrix":
    self.assertEqual(len(report["interval_ldl_pivots"]), 2)
        self.assertIn(2, report["active_prime_powers"])
        self.assertIn("not global Weil positivity", report["claim_scope"])


if __name__ == "__main__":
    unittest.main()
