"""P85 exact algebra tests; root certificate is built by separate Arb CI steps."""
import unittest
from fractions import Fraction as F
from .hilbert_source import lemma,m0_bracket,exact_boundary_cancellation
from .b4star_arb import floor_fraction,ceil_fraction,SCHEMA as ROOT_SCHEMA
from NUMERIA_P76_COEFFICIENTS.table import TABLE,TABLE_S2_PUBLIC_ROBUST

class P85Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc=lemma()

    def test_m0_rational_bounds_strict(self):
        lo,hi=map(F,m0_bracket())
        self.assertEqual(lo,F(263,200))
        self.assertEqual(hi,F(329,250))
        self.assertLess(lo,hi)

    def test_skew_cancel_exact_many(self):
        self.assertEqual(len(self.doc["finite_symbolic_pair_cancellation_tests"]),4)
        self.assertTrue(all(x["initial_boundary_exact"]=="0"
                            and x["paired_contribution_exact"]=="0"
                            for x in self.doc["finite_symbolic_pair_cancellation_tests"]))

    def test_real_coefficients_needed(self):
        with self.assertRaises(ValueError):
            exact_boundary_cancellation([1,1],[3,4])
        with self.assertRaises(ValueError):
            exact_boundary_cancellation([1,2,3],[3,4])

    def test_preissmann_not_promoted(self):
        d=self.doc
        self.assertFalse(d["source_theorem_Preissmann_proved_in_P85"])
        self.assertFalse(d["source_theorem_RZ_proved_in_P85"])
        self.assertFalse(d["X_3e10_to_3e12_covered"])
        self.assertFalse(d["global_RH_proved"])
        self.assertTrue(d["infinite_series_extension_proved_by_standard_comparison"])
        self.assertIn("pi*m0",d["elementary_transfer"]["FX_bound"])
        self.assertFalse(d["weighted_Hilbert_input"]
                         ["Carneiro_Littmann_2pi_variant_sufficient_to_keep_m0"])

    def test_floor_ceil_exact_no_floats(self):
        for v,lo,hi in ((F(7,3),2,3),(F(-7,3),-3,-2),
                        (F(5,1),5,5),(F(-2,1),-2,-2)):
            self.assertEqual(floor_fraction(v),lo)
            self.assertEqual(ceil_fraction(v),hi)

    def test_s2_table_sanity(self):
        self.assertEqual(len(TABLE),16)
        self.assertEqual(len(TABLE_S2_PUBLIC_ROBUST),16)
        for i in range(10,17):
            self.assertGreater(F(TABLE_S2_PUBLIC_ROBUST[i-1].b1),
                               F(TABLE[i-1].b1))

if __name__=="__main__":
    unittest.main()
