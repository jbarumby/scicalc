import unittest
from fractions import Fraction
from p151_dense_multiscale import (
    level_step,grid_indices,count_generators,count_local_moment_basis,
    approximate_error_bound,conservative_projection_factor,
    corrected_error_bound,local_annihilator_laurent_zero,
    fixed_level_rank_and_kernel,exact_rational_toy_moment_projection, report,
    single_scale_generator_count,single_scale_moment_dimension,
    single_scale_projection_factor,single_scale_admissible_error)

class P151Tests(unittest.TestCase):
    def test_single_scale_count(self):
        for m in range(1,10):
            self.assertEqual(single_scale_generator_count(m),2*m*2**m+1)
            self.assertEqual(single_scale_moment_dimension(m),2*m*2**m-1)
    def test_single_scale_error(self):
        self.assertEqual(single_scale_projection_factor(3),1153)
        self.assertEqual(single_scale_admissible_error(Fraction(1),Fraction(2),3,3),
                         1153*approximate_error_bound(Fraction(1),Fraction(2),3))
        with self.assertRaises(ValueError): single_scale_admissible_error(Fraction(1),Fraction(1),2,3)
    def test_step(self): self.assertEqual(level_step(3),Fraction(1,8))
    def test_grid_negative_centres(self): self.assertEqual(tuple(grid_indices(1,1)),(-2,-1,0,1,2))
    def test_formula_counts(self):
        for N in range(1,9):
            self.assertEqual(count_generators(N),2*N*(2**(N+1)-2)+N)
    def test_local_basis_counts(self):
        for N in range(1,9):
            self.assertEqual(count_local_moment_basis(N),count_generators(N)-2*N)
    def test_moment_polynomial(self):
        for j in range(-100,101): self.assertTrue(local_annihilator_laurent_zero(j))
    def test_rank_at_level(self):
        for N in range(1,5):
            for m in range(1,N+1):
                r=fixed_level_rank_and_kernel(m,N)
                self.assertTrue(r['symbolic_moment_cancellation'])
                self.assertEqual(r['kernel_dimension'],r['centres']-2)
    def test_error_m(self):
        for m in range(1,8):
            d=level_step(m)
            self.assertEqual(approximate_error_bound(Fraction(1),Fraction(1),m),
                             Fraction(5,2)*d+d*d/8)
    def test_error_decreases(self):
        es=[approximate_error_bound(Fraction(2),Fraction(3),m) for m in range(1,12)]
        self.assertTrue(all(x>y for x,y in zip(es,es[1:])))
    def test_projection_factor(self): self.assertEqual(conservative_projection_factor(3),577)
    def test_projected_error(self):
        self.assertEqual(corrected_error_bound(Fraction(1),Fraction(1),1,3),
                         577*approximate_error_bound(Fraction(1),Fraction(1),1))
    def test_exact_toy_projection(self):
        v=exact_rational_toy_moment_projection()
        self.assertEqual(len(v['projections']),5)
        self.assertTrue(all(k['residues']==['0','0'] for k in v['projections']))
    def test_report_shape(self):
        r=report(4)
        self.assertEqual(r['generator_count'],244)
        self.assertEqual(len(r['levels']),4)
    def test_bad_levels(self):
        with self.assertRaises(ValueError): level_step(0)
        with self.assertRaises(ValueError): count_generators(0)
        with self.assertRaises(ValueError): grid_indices(1,2)
    def test_negative_bound(self):
        with self.assertRaises(ValueError): approximate_error_bound(Fraction(-1),Fraction(1),1)

if __name__=='__main__':unittest.main(verbosity=2)

class P151ToeplitzTests(unittest.TestCase):
    def test_toeplitz_rational(self):
        from p151_toeplitz_exact import verify_toeplitz_congruence
        for n in range(3,15):
            w={k: Fraction(1,2)**k for k in range(n)}
            r=verify_toeplitz_congruence(w,n,Fraction(3,2))
            self.assertEqual(r['checked_entries'],(n-2)**2)
    def test_identity_lag(self):
        from p151_toeplitz_exact import gram_lag_filtered
        q=Fraction(3,2); c=q+1/q
        w={k:Fraction(int(k==0)) for k in range(11)}
        self.assertEqual(gram_lag_filtered(w,0,c),c*c+2)
        self.assertEqual(gram_lag_filtered(w,1,c),-2*c)
        self.assertEqual(gram_lag_filtered(w,2,c),1)
        self.assertEqual(gram_lag_filtered(w,3,c),0)
    def test_toeplitz_invalid(self):
        from p151_toeplitz_exact import verify_toeplitz_congruence
        with self.assertRaises(ValueError):verify_toeplitz_congruence({0:Fraction(1)},2,Fraction(2))
        with self.assertRaises(ValueError):verify_toeplitz_congruence({0:Fraction(1),1:Fraction(2),2:Fraction(3)},3,Fraction(1))
