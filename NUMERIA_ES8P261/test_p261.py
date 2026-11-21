"""Regression gates: exact finite Ising zeros, Fisher, positive definiteness and LP."""
import json
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from .leyang_fisher import (
    path_ising_coefficients, exact_path_at_zero_field, quartic_circle_criterion,
    is_palindromic, determinant_exact, rational_characteristic_grid,
    all_principal_minors_nonnegative, lp_counterexample, bost_connes_finite_gate, audit
)

class P261Tests(unittest.TestCase):
    def test_two_spin_ferromagnetic(self):
        self.assertEqual(path_ising_coefficients(2,F(1,2)),(F(1),F(1),F(1)))
    def test_three_spin_ferromagnetic(self):
        self.assertEqual(path_ising_coefficients(3,F(1,2)),
                         (F(1),F(5,4),F(5,4),F(1)))
    def test_four_spin_explicit_coefficients(self):
        self.assertEqual(path_ising_coefficients(4,F(1,2)),
                         (F(1),F(3,2),F(7,4),F(3,2),F(1)))
    def test_ferro_quartic_circle_exact(self):
        p=quartic_circle_criterion(F(1,2))
        self.assertTrue(p['all_zeros_on_unit_circle_proved'])
        self.assertEqual(p['discriminant_exact'],'13/4')
        self.assertEqual(p['Q_at_minus2_exact'],'3/4')
    def test_antiferro_quartic_offcircle_exact(self):
        p=quartic_circle_criterion(F(2))
        self.assertFalse(p['all_zeros_on_unit_circle_proved'])
        self.assertEqual(p['discriminant_exact'],'40')
        self.assertFalse(p['vertex_inside_interval'])
    def test_ferro_full_coupling_interval(self):
        for a in (F(0),F(1,4),F(1,3),F(1,2),F(3,4),F(1)):
            with self.subTest(a=a):
                self.assertTrue(quartic_circle_criterion(a)['all_zeros_on_unit_circle_proved'])
    def test_fisher_zero_field_factorization(self):
        for length in range(1,9):
            for a in (F(0),F(1,2),F(1),F(2)):
                self.assertEqual(exact_path_at_zero_field(length,a),
                                 2*(1+a)**(length-1))
    def test_spin_flip_palindromy(self):
        for length in range(1,8):
            self.assertTrue(is_palindromic(path_ising_coefficients(length,F(1,2))))
    def test_characteristic_gram_exact(self):
        m=rational_characteristic_grid(4)
        self.assertEqual(m[0],(F(1),F(3,4),F(1,2),F(3,4)))
    def test_characteristic_minors(self):
        m=rational_characteristic_grid(4)
        a=all_principal_minors_nonnegative(m)
        self.assertEqual(a['1'],['1']*4)
        self.assertEqual(a['3'],['3/16']*4)
        self.assertEqual(a['4'],['0'])
        self.assertEqual(sum(map(len,a.values())),15)
    def test_exact_determinants(self):
        self.assertEqual(determinant_exact([[2,1],[1,3]]),F(5))
        self.assertEqual(determinant_exact([[1,1],[1,1]]),F(0))
        self.assertEqual(determinant_exact([[0,1],[2,3]]),F(-2))
    def test_laguerre_local_failure(self):
        e=lp_counterexample()
        self.assertEqual(e['laguerre_inequality_at_tpi_exact'],'-1/8')
        self.assertEqual(e['laguerre_inequality_at_t0_exact'],'1/4')
        self.assertTrue(e['positive_definite_but_not_Laguerre_Polya'])
    def test_not_all_positive_coefficients_force_circle(self):
        p=path_ising_coefficients(4,F(2))
        self.assertTrue(all(c>0 for c in p))
        self.assertTrue(is_palindromic(p))
        self.assertFalse(quartic_circle_criterion(F(2))['all_zeros_on_unit_circle_proved'])
    def test_bost_connes_finite_gate(self):
        d=bost_connes_finite_gate()
        self.assertEqual(d['all_zeros_real_part_exact'],'0')
        self.assertTrue(d['critical_line_from_partition_not_proved'])
    def test_model_rejects_negative_weights(self):
        with self.assertRaises(ValueError):path_ising_coefficients(3,F(-1,2))
    def test_model_rejects_invalid_size(self):
        for n in (0,13,3.0):
            with self.assertRaises(ValueError):path_ising_coefficients(n,F(1,2))
    def test_characteristic_grid_rejects_invalid_size(self):
        for n in (0,9,2.0):
            with self.assertRaises(ValueError):rational_characteristic_grid(n)
    def test_determinant_rejects_non_square(self):
        with self.assertRaises(ValueError):determinant_exact([[1,2]])
    def test_scope_flags(self):
        a=audit()
        self.assertFalse(a['independent_Xi_Laguerre_Polya_membership_proved'])
        self.assertFalse(a['Weil_operator_positive_transfer_proved'])
        self.assertFalse(a['Riemann_hypothesis_proved'])
        self.assertEqual(a['spin_flip_and_fisher_identities_checked'],20)
    def test_json_roundtrip(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/'p261.json'
            p.write_text(json.dumps(audit(),indent=2))
            self.assertTrue(json.loads(p.read_text())['ferromagnetic_4spin_exact']['all_zeros_on_unit_circle_proved'])

if __name__=='__main__':
    unittest.main()
