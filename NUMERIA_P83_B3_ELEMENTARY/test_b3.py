"""P83 regression: exact rational proof budget and finite Möbius sanity checks.

The *general proof* is the analytic argument in README, not extrapolation
from tests at X<=200; tests certify every rational constant used there.
"""
import unittest
from fractions import Fraction as Q
from .b3_proof import (
    proof, mu_sieve, coefficient_sieve, exact_finite_short_sum,
    B3, PI_LOWER, SQFREE_ERROR_CONSTANT,
)

class ElementaryB3ProofTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=proof()

    def test_exact_global_budget_strict(self):
        d=self.d
        self.assertEqual(d["schema"],"NUMERIA.P83.INDEPENDENT_B3.SHORT_RANGE.1")
        self.assertEqual(d["X0"],10**9)
        self.assertEqual(Q(d["b3"]),Q(121,200))
        self.assertGreater(Q(d["positive_rational_slack"]),0)
        self.assertEqual(Q(d["normalized_total_majorant"])+
                         Q(d["positive_rational_slack"]),B3)
        self.assertTrue(d["b3_inequality_proved_from_elementary_standard_identities"])

    def test_machin_pi_lower_without_float(self):
        lower=16*(Q(1,5)-Q(1,375))-Q(4,239)
        self.assertEqual(Q(self.d["pi_machin_rational_lower"]),lower)
        self.assertGreater(lower,PI_LOWER)
        # tan(4 atan(1/5))=120/119; subtract atan(1/239).
        self.assertEqual(
            (Q(120,119)-Q(1,239))/(1+Q(120,119)*Q(1,239)),1
        )

    def test_squarefree_error_bound_exact(self):
        # m=floor sqrt(x), m>=r-1:
        # (m+x/m)/r <=2+1/(r-1)<=2+1/29999<2.01.
        self.assertLess(Q(2)+Q(1,29999),SQFREE_ERROR_CONSTANT)
        self.assertEqual(Q(self.d["S1_error_factor_rational_upper"]),Q(65,24))
        self.assertEqual(Q(1)+Q(3,8)+Q(4,3),Q(65,24))

    def test_mobius_sieve_consistency(self):
        mob=mu_sieve(150)
        ref={1:1,2:-1,3:-1,4:0,5:-1,6:1,8:0,
             9:0,10:1,12:0,15:1,30:-1,36:0,49:0,77:1,121:0}
        for k,v in ref.items():
            with self.subTest(k=k):self.assertEqual(mob[k],v)

    def test_all_coefficients_with_2X_to_5X_bound(self):
        # This checks the elementary case distinction; it is not proof
        # at X>=1e9, which is given symbolically in the README.
        for X in range(1,102):
            a=coefficient_sieve(X)
            self.assertTrue(all(abs(a[n])<=1 for n in range(2*X+1,5*X)))
            self.assertTrue(all(a[n]==0 for n in range(2,X+1)))

    def test_exact_truncated_sums_diagnostic(self):
        for X in (2,3,5,7,10,12,25,40):
            a=exact_finite_short_sum(X)
            self.assertTrue(a["at_most_one_for_2X_to_5X"])
            self.assertGreaterEqual(Q(a["S1"]),0)
            self.assertGreaterEqual(Q(a["S2"]),0)
            self.assertEqual(Q(a["total"]),Q(a["S1"])+Q(a["S2"]))

    def test_standard_hypotheses_not_misclaimed(self):
        d=self.d
        self.assertFalse(d["ramare_short_range_lemma_required_as_external"])
        self.assertFalse(d["mv_hilbert_inequality_independently_proved"])
        self.assertFalse(d["RZ_S_epsilon_bounds_independently_proved"])
        self.assertFalse(d["RH_globally_proved"])

    def test_adversarial_inputs_fail(self):
        with self.assertRaises(ValueError):mu_sieve(0)
        with self.assertRaises(ValueError):coefficient_sieve(-2)
        with self.assertRaises(ValueError):exact_finite_short_sum(0)

if __name__=="__main__":
    unittest.main()
