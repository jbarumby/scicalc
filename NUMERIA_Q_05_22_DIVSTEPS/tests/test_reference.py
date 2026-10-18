"""Testes de referência clássica para preparar a reprodução de divsteps.

Estes testes não simulam nem provam o circuito reversível publicado.
Executar com: python -m unittest discover -s NUMERIA_Q_05_22_DIVSTEPS/tests -v
"""
import unittest


def inverse_extended_euclid(a: int, p: int) -> int:
    """Inverso modular independente via Euclides estendido; p deve ser primo."""
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


if __name__ == "__main__":
    unittest.main()
