import unittest
from fractions import Fraction as F
from pathlib import Path
import json
from audit_es8p19 import (run, fixed_models, projected_kernel,
    correction_factors, multiply_uvt, matrix_rank, qform_on_support)

class ES8P19Tests(unittest.TestCase):
    def test_all_rational_cases(self):
        results=run()
        self.assertEqual(len(results['block_correction']['cases']),8)
        self.assertTrue(all(x['exact_outer_product_equality'] for x in results['block_correction']['cases']))
    def test_rank_bound_is_sharp(self):
        for n in (4,8,10,16,32):
            A00,A0h,Ah0,Ahh,C=fixed_models(n)
            U,V=correction_factors(A00,A0h,Ah0,C)
            self.assertEqual(matrix_rank(multiply_uvt(U,V)),4)
    def test_symmetry(self):
        for n in (3,8,16):
            A00,A0h,Ah0,Ahh,C=fixed_models(n)
            P=projected_kernel(A00,A0h,Ah0,Ahh,C)
            self.assertTrue(all(P[i][j]==P[j][i] for i in range(n) for j in range(n)))
    def test_corruption_breaks_factorization(self):
        n=8
        A00,A0h,Ah0,Ahh,C=fixed_models(n)
        P=projected_kernel(A00,A0h,Ah0,Ahh,C)
        U,V=correction_factors(A00,A0h,Ah0,C)
        E=multiply_uvt(U,V)
        P[0][0]+=F(1,1000000)
        self.assertNotEqual(P[0][0]-Ahh[0][0], E[0][0])
    def test_altered_pole_sign_detected(self):
        n=8
        A00,A0h,Ah0,Ahh,C=fixed_models(n)
        P=projected_kernel(A00,A0h,Ah0,Ahh,C)
        badC=[row[:] for row in C]
        badC[0][1]=-badC[0][1]
        U,V=correction_factors(A00,A0h,Ah0,badC)
        E=multiply_uvt(U,V)
        self.assertTrue(any(P[i][j]-Ahh[i][j]!=E[i][j] for i in range(n) for j in range(n)))
    def test_small_norm_not_finite_rank(self):
        e=F(1,1000)
        for n in (3,8,32):
            self.assertEqual(matrix_rank([[e if i==j else F() for j in range(n)] for i in range(n)]),n)
    def test_no_uniform_gap(self):
        for n in (1,4,8,16,1000):
            self.assertEqual(F(1,n),min(F(1,k+1) for k in range(n)))
        self.assertLess(F(1,1000),F(1,20))
    def test_prefix_no_extension_positivity(self):
        A=lambda i,j: F(-1) if i==j==9 else F(1) if i==j else F()
        for i in range(8): self.assertEqual(qform_on_support(A,{i:F(1)}),1)
        self.assertEqual(qform_on_support(A,{9:F(1)}),-1)
    def test_support_insertion_and_equality(self):
        A=lambda i,j:F(2) if i==j else F(1,7)
        x={3:F(2),11:F(-1),200:F(1,3)}
        y={11:F(-1),200:F(1,3),3:F(2)}
        self.assertEqual(qform_on_support(A,x),qform_on_support(A,y))
        self.assertEqual(qform_on_support(A,{200:F(1,3)}),F(2,9))
    def test_bad_quoted_qwen_implication(self):
        # Rank of each finite truncation ≤n doesn't yield a uniform rank bound.
        ranks=[matrix_rank([[F(1,1000) if i==j else F() for j in range(n)] for i in range(n)]) for n in (4,8,12)]
        self.assertEqual(ranks,[4,8,12])
    def test_reproduction(self):
        p=Path(__file__).with_name('RESULTADOS_ES8P19.json')
        self.assertEqual(json.loads(p.read_text()),run())

if __name__=='__main__':unittest.main()
