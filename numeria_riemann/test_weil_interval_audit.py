
    def test_small_interval_audit_returns_explicit_status(self):
        report = run(
            L_value=0.45,
            N=4,
            R_value=40,
            step_value=0.05,
            terms=32,
            precision=128,
        )

        self.assertIn(report["status"], {
            "certified_positive_definite_finite_matrix",
            "inconclusive_interval_ldl",
        })

        self.assertEqual(report["projected_dimension"], 2)
        self.assertEqual(report["parameters"]["quadrature_cells"], 800)

        pivots = report["interval_ldl_pivots"]
        self.assertGreaterEqual(len(pivots), 1)
        self.assertLessEqual(len(pivots), 2)

        if report["status"] == "certified_positive_definite_finite_matrix":
            self.assertEqual(len(pivots), 2)

        self.assertIn(2, report["active_prime_powers"])
        self.assertIn(
            "not global Weil positivity",
            report["claim_scope"],
        )
