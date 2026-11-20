"""Gate-level block encoding of the P249 prime shifts on uniform cells."""

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from math import fsum
from operator import index
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, qasm2, transpile
from qiskit import __version__ as qiskit_version
from qiskit.circuit.library import QFTGate, StatePreparation, U3Gate, UCGate
from qiskit.quantum_info import Statevector

from NUMERIA_ES8P249.prime_operator import PrimeShiftOperator


HALF_WIDTH = Fraction(19, 20)


@dataclass(frozen=True)
class PrimeLCU:
    cells: int
    padded_size: int
    coefficients: tuple[float, ...]
    shifts: tuple[int, ...]
    normalization: float
    data_qubits: int
    selector_qubits: int


def _flip(circuit, controls, target):
    if controls:
        circuit.mcx(controls, target)
    else:
        circuit.x(target)


def _select_shift(circuit, data, selector, branch, shift):
    """Apply a constant modular translation on one selector basis state."""
    if shift == 0:
        return
    zero_bits = [qubit for bit, qubit in enumerate(selector) if branch & (1 << bit) == 0]
    for qubit in zero_bits:
        circuit.x(qubit)

    # Bitwise complement conjugates addition to subtraction modulo 2**len(data).
    if shift < 0:
        for qubit in data:
            circuit.x(qubit)
    for bit in range(len(data)):
        if abs(shift) & (1 << bit):
            # Increment the upper register; controls see its original low bits.
            for target in range(len(data) - 1, bit, -1):
                _flip(circuit, selector + data[bit:target], data[target])
            _flip(circuit, selector, data[bit])
    if shift < 0:
        for qubit in data:
            circuit.x(qubit)

    for qubit in zero_bits:
        circuit.x(qubit)


def _select_fourier(circuit, data, selector, plan):
    if not any(plan.shifts):
        return
    padded_shifts = list(plan.shifts) + [0] * ((1 << len(selector)) - len(plan.shifts))
    fourier = QFTGate(len(data))
    circuit.append(fourier, data)
    for bit, qubit in enumerate(data):
        angles = [2 * np.pi * ((shift << bit) % plan.padded_size) / plan.padded_size
                  for shift in padded_shifts]
        if selector:
            gates = [np.diag([1, np.exp(1j * angle)]) for angle in angles]
            circuit.append(UCGate(gates), [qubit, *selector])
        else:
            circuit.p(angles[0], qubit)
    circuit.append(fourier.inverse(), data)


def build_prime_lcu(cells, *, half_width=HALF_WIDTH, select_method="fourier"):
    """Return PREP–SELECT–UNPREP for the four P249 prime powers.

    The first cells basis states form the physical register subspace. The
    ancilla-zero block compressed to those states is the Toeplitz operator
    divided by plan.normalization. The remaining data states implement
    padding and are part of the postselection condition.
    """
    if select_method not in ("fourier", "reversible"):
        raise ValueError("select_method must be fourier or reversible")
    cells = index(cells)
    if cells < 1:
        raise ValueError("cells must be positive")
    lags = PrimeShiftOperator(cells, half_width).lags
    terms = [(0, float(lags[0]))] if lags[0] > 0 else []
    for shift, weight in enumerate(lags[1:], start=1):
        if weight > 0:
            terms.extend(((shift, float(weight)), (-shift, float(weight))))
    if not terms:
        raise ValueError("the selected prime shift operator is zero")

    shifts, coefficients = zip(*terms, strict=True)
    normalization = fsum(coefficients)
    data_qubits = max(1, (2 * cells - 2).bit_length())
    selector_qubits = (len(terms) - 1).bit_length()
    plan = PrimeLCU(
        cells, 1 << data_qubits, coefficients, shifts, normalization, data_qubits, selector_qubits,
    )
    circuit = QuantumCircuit(data_qubits + selector_qubits, name="prime_toeplitz_lcu")
    data = list(range(data_qubits))
    selector = list(range(data_qubits, circuit.num_qubits))
    if selector:
        state = np.zeros(1 << selector_qubits)
        state[:len(terms)] = np.sqrt(np.array(coefficients) / normalization)
        prepare = StatePreparation(state, normalize=True)
        circuit.append(prepare, selector)

    if select_method == "fourier":
        _select_fourier(circuit, data, selector, plan)
    else:
        for branch, shift in enumerate(shifts):
            _select_shift(circuit, data, selector, branch, shift)
    if selector:
        circuit.append(prepare.inverse(), selector)
    return circuit, plan


def _portable_compile(circuit):
    compiled = transpile(circuit, basis_gates=["u3", "cx"], optimization_level=1, seed_transpiler=17)
    # OpenQASM 2 has no global phase declaration. Encode P(phi) X P(phi) X
    # explicitly so the exported ancilla-zero block retains its phase.
    phase = float(compiled.global_phase)
    if phase:
        compiled.global_phase = 0
        for _ in range(2):
            compiled.append(U3Gate(0, 0, phase), [0])
            compiled.append(U3Gate(np.pi, 0, np.pi), [0])
    return compiled


def simulate_demo(cells, *, half_width=HALF_WIDTH, select_method="fourier"):
    """Run a deterministic complex input and compile to a portable gate set."""
    circuit, plan = build_prime_lcu(cells, half_width=half_width, select_method=select_method)
    indices = np.arange(cells)
    vector = ((5 * indices + 1) % 13 - 6) + 1j * ((7 * indices + 2) % 11 - 5)
    vector = vector / np.linalg.norm(vector)
    initial = np.zeros(1 << circuit.num_qubits, dtype=complex)
    initial[:cells] = vector
    output = Statevector(initial).evolve(circuit).data
    selected = output[:cells]
    expected = PrimeShiftOperator(cells, half_width).matvec(vector)
    padding = output[cells:plan.padded_size]
    compiled = _portable_compile(circuit)
    exported = qasm2.loads(qasm2.dumps(compiled))
    archive_output = Statevector(initial).evolve(exported).data
    return {
        "schema": "NUMERIA.ES8P252.PRIME_SHIFT_LCU.1",
        "qiskit_version": qiskit_version,
        "numpy_version": np.__version__,
        "half_width_exact": str(Fraction(half_width)),
        "cells": cells,
        "select_method": select_method,
        "padded_size": plan.padded_size,
        "data_qubits": plan.data_qubits,
        "selector_qubits": plan.selector_qubits,
        "shifts": plan.shifts,
        "coefficients_binary64": plan.coefficients,
        "LCU_normalization_binary64": plan.normalization,
        "statevector_max_abs_error_binary64": float(np.max(abs(plan.normalization * selected - expected))),
        "OpenQASM2_round_trip_max_abs_error_binary64": float(np.max(abs(archive_output - output))),
        "probability_selector_zero_and_physical_cells": float(np.vdot(selected, selected).real),
        "probability_selector_zero_and_padding": float(np.vdot(padding, padding).real),
        "compiled_gate_counts": {name: int(count) for name, count in compiled.count_ops().items()},
        "compiled_depth": compiled.depth(),
    }, compiled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=16)
    parser.add_argument("--half-width", type=Fraction, default=HALF_WIDTH)
    parser.add_argument("--select", choices=("fourier", "reversible"), default="fourier")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qasm", type=Path)
    args = parser.parse_args()
    report, circuit = simulate_demo(args.cells, half_width=args.half_width, select_method=args.select)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if args.qasm is not None:
        args.qasm.parent.mkdir(parents=True, exist_ok=True)
        args.qasm.write_text(qasm2.dumps(circuit) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
