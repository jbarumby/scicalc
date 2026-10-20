"""Independent classical reference for half-delta divsteps (not a circuit proof)."""
import random
import unittest

P = 2**256 - 2**32 - 977
BOUND = 590
CIRCUIT_SCHEDULE = 533
INV2 = (P + 1) // 2


def half(v: int) -> int:
    return (v * INV2) % P


def hddivsteps(x: int, max_steps: int = BOUND + 1) -> tuple[int, int]:
    if not 0 <= x < P:
        raise ValueError("x must be canonical: 0 <= x < P")
    delta, f, g, d, e = 1, P, x, 0, 1
    steps = 0
    while g != 0:
        if steps >= max_steps:
            raise AssertionError("walk exceeded max_steps")
        if g & 1 and delta > 0:
            delta, f, g, d, e = 2 - delta, g, (g - f) // 2, e, half(e - d)
        elif g & 1:
            delta, f, g, d, e = 2 + delta, f, (g + f) // 2, d, half(e + d)
        else:
            delta, f, g, d, e = 2 + delta, f, g // 2, d, half(e)
        steps += 1
    if x == 0:
        return steps, 0
    if abs(f) != 1:
        raise AssertionError(f"terminal gcd invariant failed: f={f}")
    return steps, (d * f) % P


def edge_inputs() -> list[int]:
    values = {0, 1, 2, 3, P - 1, P - 2, P - 3, (P + 1) // 2, (P - 1) // 2}
    values.update(2**k for k in (1, 31, 32, 64, 128, 255))
    values.update(2**k - 1 for k in (2, 32, 64, 128, 255))
    values.update((2**256 - 2**32 - 978, 2**32, 2**32 + 977))
    return sorted(x for x in values if 0 <= x < P)


class TestHalfDeltaDivstepsReference(unittest.TestCase):
    def test_modulus_is_secp256k1_prime(self):
        self.assertEqual(P, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F)

    def test_edge_inputs(self):
        for x in edge_inputs():
            steps, inv = hddivsteps(x)
            if x == 0:
                self.assertEqual((steps, inv), (0, 0))
            else:
                self.assertEqual(inv, pow(x, -1, P), x)
                self.assertLessEqual(steps, BOUND, x)

    def test_random_inputs_match_builtin_inverse_and_bound(self):
        rng = random.Random(20261018)
        n = 1500
        for _ in range(n):
            x = rng.randrange(1, P)
            steps, inv = hddivsteps(x)
            self.assertEqual((x * inv) % P, 1)
            self.assertLessEqual(steps, BOUND)

    def test_rejects_noncanonical_inputs(self):
        for x in (-1, P, P + 1):
            with self.assertRaises(ValueError):
                hddivsteps(x)

    def test_inverse_property_not_a_circuit_test(self):
        _, inv = hddivsteps(12345)
        self.assertEqual((12345 * inv) % P, 1)
        self.assertNotEqual((12346 * inv) % P, 1)


if __name__ == "__main__":
    unittest.main()
