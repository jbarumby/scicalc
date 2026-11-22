"""P92G arithmetic regression, all exact Fraction, no assertions as proof."""
import unittest
from fractions import Fraction as F
from .certificate import (
    SCHEMA, TARGET, certify, points_and_deltas, entries,
    geometric_kernel, kernel_one_side_ceiling, audit_case,
    exact_ldlt_positive
)


class P92GTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=certify()

    def test_reproducibility_and_scientific_scope(self):
        d=self.data
        self.assertEqual(d["schema"],SCHEMA)
        self.assertEqual(d["rational_regressions"],57)
        self.assertEqual(d["independent_rational_ldlt"],42)
        self.assertTrue(d["all_real_geometric_ratios_proved_analytically"])
        self.assertTrue(d["number_of_distinct_scales_unbounded_in_theorem"])
        for field in ("arbitrary_point_configurations_proved",
                      "general_Preissmann_Lemma6_proved",
                      "Riemann_Hypothesis_proved"):
            self.assertFalse(d[field],field)

    def test_geometric_distances(self):
        for r in (F(101,100),F(6,5),F(3,2),F(2),F(7)):
            for n in (2,3,8,15):
                x,d=points_and_deltas(r,n)
                self.assertEqual(x[0],0)
                self.assertEqual(x[1],1)
                self.assertEqual([x[k+1]-x[k] for k in range(n-1)],
                                 [r**(2*k) for k in range(n-1)])
                self.assertEqual(d,[F(1)]+[r**(2*(k-1)) for k in range(1,n)])
                self.assertEqual(d,[min(abs(x[i]-x[j]) for j in range(n) if j!=i)
                                    for i in range(n)])

    def test_exact_toeplitz_rows(self):
        for r in (F(11,10),F(3,2),F(2),F(3)):
            M=entries(r,12)
            for i in range(1,12):
                for j in range(1,12):
                    if i==j:
                        self.assertEqual(M[i][j],0)
                    else:
                        self.assertEqual(M[i][j],geometric_kernel(r,abs(i-j)))
                if i>=2:
                    self.assertGreaterEqual(M[1][i],M[0][i])
            self.assertEqual(M[0][1],0)

    def test_kernel_series_ceiling(self):
        for r in (F(1001,1000),F(11,10),F(3,2),F(2),F(5)):
            S=kernel_one_side_ceiling(r)
            for kmax in (2,10,30):
                finite=sum((geometric_kernel(r,k) for k in range(1,kmax)),F(0))
                self.assertLessEqual(finite,S)
            self.assertLess(3*S,TARGET)

    def test_quartic_maximum_has_strict_rational_ceiling(self):
        self.assertGreater(F(71,17)**2,17)
        self.assertEqual((F(107)+51*F(71,17))/512,F(5,8))
        for n in (2,3,4,7,10,25,100):
            for k in range(n+1):
                t=F(k,n)
                f=t*(1-t)*(1+t)**2
                self.assertLessEqual(f,F(5,8))
                self.assertEqual(f,t+t*t-t**3-t**4)

    def test_all_exact_regressions_below_threshold(self):
        for item in self.data["cases"]:
            val=F(item["actual_max_row_sum"])
            bound=F(item["universal_rational_ceiling"])
            self.assertLessEqual(val,bound)
            self.assertLess(bound,TARGET)
            self.assertEqual(F(item["margin_to_15_over_16"]),TARGET-val)

    def test_independent_ldlt(self):
        for r in (F(6,5),F(2),F(3)):
            M=entries(r,8)
            pivot=exact_ldlt_positive(TARGET,M)
            self.assertGreater(pivot,0)

    def test_invalid_parameters(self):
        for r,n in ((1,3),(F(1,2),4),(F(3,2),1)):
            with self.assertRaises(ValueError):
                points_and_deltas(r,n)
        with self.assertRaises(ValueError):
            geometric_kernel(F(3,2),0)
        with self.assertRaises(ValueError):
            kernel_one_side_ceiling(1)


if __name__=="__main__":
    unittest.main()
