"""Testes de referência clássica para preparar a reprodução de divsteps.

Estes testes não simulam nem provam o circuito reversível publicado.
Executar com: python -m unittest discover -s NUMERIA_Q_05_22_DIVSTEPS/tests -v
"""
import unittest


def inverse_extended_euclid(a: int, p: int) -> int:
    """Inverso modular independente via Euclides estendido; gcd(a, p) deve ser 1."""
    old_r, r = a, p
    old_s, s = 1, 0
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    if old_r != 1:
        raise ValueError("inverso modular inexistente")
    return old_s % p


class TestClassicalReference(unittest.TestCase):
    def test_all_nonzero_residues_small_primes(self):
        for p in (2, 3, 5, 7, 11, 13, 17, 19):
            for a in range(1, p):
                inv = inverse_extended_euclid(a, p)
                self.assertEqual((a * inv) % p, 1, (a, p, inv))

    def test_zero_has_no_inverse(self):
        for p in (2, 3, 5, 7, 11):
            with self.assertRaises(ValueError):
                inverse_extended_euclid(0, p)

    def test_out_of_range_non_coprime_inputs_rejected(self):
        for a, p in ((2, 4), (6, 9), (0, 13)):
            with self.assertRaises(ValueError):
                inverse_extended_euclid(a, p)

    def test_large_modular_inverse_edge_case(self):
        # Para a = 2**59 - 1 e p = 2**59 + 1, a = p - 2.
        # O inverso é 2**58, pois (p - 2) * 2**58 = 2**59 * 2**58
        # é congruente a 1 módulo p, já que 2**59 = p - 1.
        a = 2**59 - 1
        p = 2**59 + 1
        expected_result = 2**58
        self.assertEqual(inverse_extended_euclid(a, p), expected_result)
        self.assertEqual((a * expected_result) % p, 1)

    def test_example_with_negative_intermediate_coefficients(self):
        # Euclides estendido pode produzir coeficientes intermediários negativos.
        # Para 3 modulo 7, o inverso é 5. Isso não inspeciona o delta interno
        # do algoritmo divsteps; testa apenas a referência modular clássica.
        a, p = 3, 7
        expected_result = 5
        self.assertEqual(inverse_extended_euclid(a, p), expected_result)
        self.assertEqual((a * expected_result) % p, 1)


if __name__ == "__main__":
    unittest.main()
