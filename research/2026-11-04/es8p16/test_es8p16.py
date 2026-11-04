import unittest
from pathlib import Path
from fractions import Fraction
from tempfile import TemporaryDirectory
from audit_es8p16 import parse_matrix, verify, DEFAULT_SOURCE

class TestES8P16(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=DEFAULT_SOURCE.read_text(encoding='utf-8')
        cls.results=verify(DEFAULT_SOURCE)

    def test_all_eight_nestings(self):
        self.assertEqual(len(self.results['nested_exact_rational_comparisons']),8)
        self.assertEqual(sum(self.results['nested_exact_rational_comparisons'].values()),208)

    def test_frobenius_finite_not_infinite(self):
        for n in [4,6,8]:
            d=self.results['frobenius_exploratory_finite'][str(n)]
            self.assertLess(Fraction(d['frobenius_norm_squared']),Fraction(1,10000))
        self.assertNotIn('infinite_frobenius_bound',self.results)

    def test_exact_row_budgets(self):
        expected=[Fraction(1937,10**6),Fraction(4379,10**6),Fraction(8333,10**6)]
        self.assertEqual([Fraction(self.results['error_row_budgets'][str(n)]['max_row_fraction']) for n in [4,6,8]],expected)

    def test_counterexample_spectral(self):
        c=self.results['ldl_pivots_not_spectral_gap']
        self.assertEqual(Fraction(c['rayleigh_fraction']),Fraction(1,101))
        self.assertLess(Fraction(c['rayleigh_fraction']),Fraction(1,20))

    def test_counterexample_projection(self):
        self.assertEqual(Fraction(2)**2*Fraction(1,100),Fraction(1,25))

    def test_shift_pivots_not_additive(self):
        p=self.results['shifted_ldl_pivot_nonadditivity']
        self.assertNotEqual(Fraction(p['pivot_at_index2_of_center'])-Fraction(p['pivot_at_index2_of_shifted']),Fraction(3,50))

    def test_symbolic_two_extensions(self):
        self.assertTrue(self.results['two_distinct_infinite_extensions_agree_on_8x8'])

    def test_mutation_detected_of_nested_error(self):
        mutated=self.source.replace('def error6 (i j : Fin 6) : ℚ :=\n  match i.val, j.val with\n  | 0, 0 => (139 : ℚ) / 500000',
                'def error6 (i j : Fin 6) : ℚ :=\n  match i.val, j.val with\n  | 0, 0 => (140 : ℚ) / 500000')
        self.assertNotEqual(mutated,self.source)
        with TemporaryDirectory() as td:
            p=Path(td)/'corrupt.lean';p.write_text(mutated,encoding='utf-8')
            with self.assertRaises(AssertionError): verify(p)

    def test_mutation_detected_missing_matrix(self):
        mutated=self.source.replace('def center8 (i j : Fin 8)', 'def destroyed8 (i j : Fin 8)')
        with TemporaryDirectory() as td:
            p=Path(td)/'corrupt.lean';p.write_text(mutated,encoding='utf-8')
            with self.assertRaises(ValueError): verify(p)

if __name__=='__main__': unittest.main()
