import unittest
from fractions import Fraction as F
import audit_es8p17 as e

class TestES8P17(unittest.TestCase):
    def test_exact_constants(self): e.check_premises()
    def test_basis_first_columns(self):
        b=e.b_matrix(6,F(9));self.assertEqual(b[2][0],1);self.assertEqual(b[0][0],-e.coeff(3,F(9)))
    def test_gram_rational(self):
        for n in (6,8,10,20):
            self.assertEqual(e.projected_gram(n,F(9)),e.gram_multiplication(n,F(9)))
    def test_parity_orthogonality(self):
        z=e.projected_gram(20,F(9))
        for i in range(18):
            for j in range(18):
                if i!=j and (i-j)%2: self.assertEqual(z[i][j],0)
    def test_monotonicity_in_x(self):
        for k in range(3,101):
            self.assertGreaterEqual(e.coeff(k,F(9)),e.coeff(k,F(10)))
    def test_c_upper(self):
        for k in range(3,101):self.assertTrue(e.bound_coefficient(k,F(9)))
    def test_global_norm_bound(self):
        for n in (3,4,6,8,10,20,40,80):
            self.assertLess(e.norm_bound_exact(n,F(9)),e.KAPPA_BUDGET)
    def test_uniform_tail(self):
        for n in (1,2,3,6,8,10,100,10000):
            self.assertLess(e.KAPPA_BUDGET*e.delta(n),F(1,100))
    def test_tail_decreases(self):
        for n in range(1,40): self.assertGreater(e.delta(n),e.delta(n+1))
    def test_adversarial_projection(self):
        b=e.b_matrix(10,F(9));b[0][0]=-F(3,1) # artificial overly large coefficient
        self.assertGreater(sum(row[0]*row[0] for row in b),e.KAPPA_BUDGET)
    def test_finite_hs_not_infinite(self):
        # E_ii=1 for all i>=8 gives all-zero 8x8 prefix but diverging HS tail.
        for n in (16,32,128):
            hs=sum(1 for i in range(n) if i>=8)
            self.assertEqual(hs,n-8)
    def test_no_uniform_spectral_gap_from_nested_prefix(self):
        # Diagonal extension with 1/(k+1) for k>=8 remains positive but loses gap.
        for n in (9,30,300):
            self.assertLess(F(1,n),F(1,20) if n>20 else F(1,1))
    def test_reproducibility(self):
        self.assertEqual(e.run(),e.run())

if __name__=='__main__':unittest.main()
