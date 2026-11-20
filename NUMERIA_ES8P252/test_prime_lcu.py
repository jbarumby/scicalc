"""Statevector regressions for the prime-only Toeplitz block encoding."""

from fractions import Fraction
import math
import unittest

import numpy as np
from qiskit import qasm2
from qiskit.quantum_info import Operator, Statevector

from NUMERIA_ES8P252.prime_lcu import _portable_compile, build_prime_lcu


def dense_overlap_reference(cells, half_width):
    """Integrate products of normalized indicators by endpoint intersections."""
    length = 2 * float(half_width)
    edges = np.linspace(-length / 2, length / 2, cells + 1)
    cell_width = length / cells
    prime_powers = ((2, 2), (3, 3), (4, 2), (5, 5))
    result = np.zeros((cells, cells))
    for row in range(cells):
        for column in range(cells):
            for power, prime in prime_powers:
                coefficient = math.log(prime) / math.sqrt(power)
                for shift in (-math.log(power), math.log(power)):
                    start = max(edges[row], edges[column] + shift)
                    stop = min(edges[row + 1], edges[column + 1] + shift)
                    result[row, column] += coefficient * max(0, stop - start) / cell_width
    return result


def _evolve_supported(circuit, vector):
    state = np.zeros(1 << circuit.num_qubits, dtype=complex)
    state[:len(vector)] = vector
    return Statevector(state).evolve(circuit).data


class ToeplitzLCUTests(unittest.TestCase):
    def test_full_selector_zero_block_matches_cell_overlap_integrals(self):
        half_width = Fraction(19, 20)
        cells = 4
        circuit, plan = build_prime_lcu(cells, half_width=half_width)
        self.assertEqual(plan.cells, cells)
        self.assertEqual(plan.padded_size, 1 << plan.data_qubits)
        self.assertGreaterEqual(plan.padded_size, 2 * cells - 1)
        self.assertEqual(circuit.num_qubits, plan.selector_qubits + plan.data_qubits)
        self.assertEqual(circuit.num_clbits, 0)
        self.assertFalse({"measure", "reset"}.intersection(circuit.count_ops()))
        self.assertEqual(len(plan.coefficients), len(plan.shifts))
        self.assertLessEqual(len(plan.shifts), 1 << plan.selector_qubits)
        self.assertTrue(all(coefficient > 0 for coefficient in plan.coefficients))
        self.assertAlmostEqual(float(plan.normalization), math.fsum(plan.coefficients), places=13)

        operator = Operator(circuit)
        self.assertTrue(operator.is_unitary())
        block = float(plan.normalization) * operator.data[:cells, :cells]
        np.testing.assert_allclose(block, dense_overlap_reference(cells, half_width), atol=2e-13, rtol=3e-13)

    def test_complex_input_action_and_postselection_probability(self):
        half_width = Fraction(19, 20)
        vector = np.array([1 + 2j, -0.2 + 1j, 0.3j, 4 - 1j, -0.4])
        vector /= np.linalg.norm(vector)
        circuit, plan = build_prime_lcu(len(vector))
        evolved = _evolve_supported(circuit, vector)
        projected = evolved[:len(vector)]
        expected = dense_overlap_reference(len(vector), half_width) @ vector
        np.testing.assert_allclose(float(plan.normalization) * projected, expected, atol=2e-13, rtol=3e-13)
        success = np.vdot(projected, projected).real
        self.assertAlmostEqual(success, np.vdot(expected, expected).real / float(plan.normalization)**2, places=13)

    def test_endpoint_impulse_requires_projection_without_circular_wrap(self):
        half_width = Fraction(19, 20)
        cells = 5
        circuit, plan = build_prime_lcu(cells)
        vector = np.zeros(cells)
        vector[-1] = 1
        projected = _evolve_supported(circuit, vector)
        expected = dense_overlap_reference(cells, half_width)[:, -1]
        np.testing.assert_allclose(float(plan.normalization) * projected[:cells], expected, atol=2e-13, rtol=3e-13)
        self.assertGreater(np.linalg.norm(projected[cells:plan.padded_size]), 0.05)
        self.assertGreater(abs(expected[0] - expected[-2]), 0.02)

    def test_one_cell_encodes_constant_overlap(self):
        circuit, plan = build_prime_lcu(1, half_width=Fraction(2, 5))
        self.assertEqual(plan.padded_size, 2)
        self.assertEqual(plan.selector_qubits, 0)
        scalar = _evolve_supported(circuit, np.ones(1))[0]
        expected = 2 * math.log(2) / math.sqrt(2) * (1 - math.log(2) / 0.8)
        self.assertAlmostEqual((float(plan.normalization) * scalar).real, expected, places=13)
        self.assertAlmostEqual(scalar.imag, 0, places=13)

    def test_invalid_or_zero_overlap_parameters_fail_closed(self):
        with self.assertRaises(ValueError):
            build_prime_lcu(0)
        with self.assertRaises(ValueError):
            build_prime_lcu(5, half_width=Fraction(1, 4))

    def test_both_select_implementations_and_qasm_preserve_the_complex_phase(self):
        vector = np.array([1j, -1, 0.5 + 0.25j])
        vector /= np.linalg.norm(vector)
        fourier, _ = build_prime_lcu(len(vector), select_method="fourier")
        reversible, _ = build_prime_lcu(len(vector), select_method="reversible")
        expected = _evolve_supported(fourier, vector)
        np.testing.assert_allclose(_evolve_supported(reversible, vector), expected, atol=1e-12)
        for circuit in (fourier, reversible):
            restored = qasm2.loads(qasm2.dumps(_portable_compile(circuit)))
            np.testing.assert_allclose(_evolve_supported(restored, vector), expected, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
