"""P87 fail-closed unit tests: Shan convexity vs Lemma6 finite witnesses."""
import unittest
from fractions import Fraction as F

from .shan_spacing import (
    F_N, mass_bound, lift_certificate, harmonic2,
    audit_nodes,proof_ledger,SCHEMA as SHAN_SCHEMA,
)
from .lemma6_bridge import (
    atan_partial, machin_pi_lower, lemma6_finite,quadratic_surd_chain,
    audit as lemma6_audit,SCHEMA as L6_SCHEMA
)
from NUMERIA_P86_PREISSMANN.spectral_rational import row_energies

class P87Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.shan=proof_ledger()
        cls.l6=lemma6_audit()

    def test_all_nodes_and_exact_normalization(self):
        self.assertEqual(self.shan["schema"],SHAN_SCHEMA)
        self.assertEqual(self.shan["exact_case_count"],9)
        self.assertTrue(self.shan["Preissmann_Lemma1_reconstructed"])
        for cert in self.shan["cases"]:
            nodes=cert["nodes"]
            energies=row_energies(nodes)
            self.assertEqual([str(v) for v in energies],
                             [r["row_energy"] for r in cert["rows"]])
            for row in cert["rows"]:
                for side in row["sides"]:
                    if side.get("empty_side"):continue
                    witness=side["gap_lift"]
                    self.assertLessEqual(F(witness["start"]),
                                         F(witness["after_lift"]))
                    self.assertLessEqual(F(witness["after_lift"]),
                                         F(witness["mass"]["mass_left"]))
                    self.assertLessEqual(F(witness["mass"]["mass_left"]),
                                         F(witness["mass"]["finite_lattice_budget"]))

    def test_adversarial_gap_descents_are_lifted(self):
        for gaps in (
            ["1","1/50","4","1/100","3"],
            ["1","5","1/11","7/2","1/17"],
            ["100","1/100","1/1000","1/2"],
        ):
            certificate=lift_certificate(gaps)
            self.assertTrue(certificate["all_exact_rational_gates_pass"])
            self.assertGreaterEqual(len(certificate["descent_steps"]),1)
            self.assertTrue(all(s["monotone_step"]
                                for s in certificate["descent_steps"]))

    def test_convex_lattice_budget_with_fractional_endpoint(self):
        val=mass_bound([F(7,5),F(17,10),F(3),F(3)])
        self.assertEqual(val["floor_end"],6)
        self.assertEqual(F(val["fractional_end"]),F(1,10))
        self.assertLess(F(val["mass_left"]),F(val["finite_lattice_budget"]))

    def test_Basel_partial_and_uniform_grid(self):
        for n in range(1,31):
            self.assertLess(harmonic2(n),F(2)) # Basel = π²/6 <2
        sample=audit_nodes(range(17))
        self.assertEqual(len(sample["rows"]),17)
        self.assertEqual(F(sample["rows"][8]["row_energy"]),
                         2*harmonic2(8))

    def test_reject_zero_degenerate_or_nonpositive(self):
        for a in ((0,1),(1,0),(F(1,2),F(1)),(1,)):
            with self.subTest(gaps=a),self.assertRaises(ValueError):
                lift_certificate(a)
        with self.assertRaises(ValueError):
            audit_nodes([1,1,2])
        with self.assertRaises(ValueError):
            mass_bound([1,F(1,2),1])

    def test_lemma6_samples_are_finite_only(self):
        self.assertEqual(self.l6["schema"],L6_SCHEMA)
        self.assertGreater(self.l6["finite_exact_pair_count"],100)
        self.assertFalse(self.l6["Preissmann_Lemma6_all_N_proved_here"])
        for witness in self.l6["finite_samples"].values():
            self.assertGreater(F(witness["minimum_exact_slack"]),0)
            self.assertTrue(witness["all_finite_gates_pass"])
            self.assertEqual(len(witness["pairs"]),witness["number_pairs"])
            for row in witness["pairs"]:
                self.assertEqual(F(row["rhs_with_pi_lower_exact"])-F(row["left_exact"]),
                                 F(row["slack_exact"]))
                self.assertGreater(F(row["slack_exact"]),0)

    def test_machin_pi_lower_exact(self):
        lo=machin_pi_lower()
        self.assertGreater(lo,F(157,50))
        self.assertLess(lo,F(22,7))
        self.assertEqual(atan_partial(F(1,5),1),F(1,5)-F(1,375))

    def test_conditional_m0_polynomial_identity(self):
        q=quadratic_surd_chain()
        self.assertTrue(q["algebra_closed_given_Lemma6"])
        self.assertEqual(q["exact_quadratic_identity_lhs"],["1/45","0"])
        self.assertEqual(q["exact_m0_squared"],["1","2/3"])

    def test_no_false_full_proof_or_repeated_S3(self):
        self.assertFalse(self.shan["Preissmann_Lemma6_reconstructed"])
        self.assertFalse(self.shan["Preissmann_full_m0_reproved_without_Lemma6"])
        self.assertFalse(self.l6["Preissmann_constant_m0_unconditionally_independent"])
        self.assertFalse(self.l6["C3_C4_B1_recertification_repeated"])
        self.assertFalse(self.shan["RH_proved"])
        self.assertFalse(self.l6["global_RH_proved"])

if __name__=="__main__":
    unittest.main()
