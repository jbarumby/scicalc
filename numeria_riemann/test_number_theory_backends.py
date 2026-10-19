import math
import unittest

from numeria_riemann.number_theory_backends import sympy_von_mangoldt


class NumberTheoryBackendTests(unittest.TestCase):
    def test_sympy_von_mangoldt_when_available(self):
        try:
            import sympy  # noqa: F401
        except ImportError:
            self.skipTest("optional SymPy dependency not installed")
        self.assertAlmostEqual(sympy_von_mangoldt(8), math.log(2))
        self.assertAlmostEqual(sympy_von_mangoldt(9), math.log(3))
        self.assertEqual(sympy_von_mangoldt(12), 0.0)
        self.assertEqual(sympy_von_mangoldt(1), 0.0)


if __name__ == "__main__":
    unittest.main()
