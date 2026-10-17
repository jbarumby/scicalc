#!/usr/bin/env python3
"""Avaliador independente: equivalência exata e recursos lógicos."""
from __future__ import annotations
import hashlib
from numq_model import Circuit, run_circuit
from classical_reference import modular_inverse_euclid

def reference_function(x: int, p: int, n: int) -> int:
    if 1 <= x < p: return modular_inverse_euclid(x, p)
    return 0

def semantic_hash(p: int, n: int) -> str:
    payload = bytes(reference_function(x, p, n) for x in range(1 << n))
    return hashlib.sha256(f"modinv-prime-p={p};n={n};f(0)=f(x>=p)=0;".encode() + payload).hexdigest()

def verify_exact(circuit: Circuit, p: int) -> dict:
    n = circuit.n_input
    if circuit.n_output != n:
        return {"passed": False, "reason": "benchmark exige registradores x/y de mesmo tamanho"}
    mask = (1 << n) - 1; checks = 0
    for x in range(1 << n):
        f = reference_function(x, p, n)
        for y in range(1 << n):
            state = x | (y << n)
            final = run_circuit(state, circuit)
            x_out = final & mask; y_out = (final >> n) & mask; anc_out = final >> (2*n)
            checks += 1
            if x_out != x or y_out != (y ^ f) or anc_out != 0:
                return {"passed": False, "reason": "mismatch", "x": x, "y": y,
                        "expected_x": x, "expected_y": y ^ f, "actual_x": x_out,
                        "actual_y": y_out, "ancilla_out": anc_out, "checks_before_failure": checks}
    return {"passed": True, "checks": checks, "inputs_x": 1 << n,
            "targets_y": 1 << n, "ancillas_clean": True,
            "specification": "|x>|y>|0> -> |x>|y XOR f(x)>|0>"}

def resource_metrics(circuit: Circuit) -> dict:
    counts = {"X": 0, "CX": 0, "CCX": 0}; last_layer = [0] * circuit.n_qubits
    for gate in circuit.gates:
        counts[gate.name] += 1
        layer = 1 + max(last_layer[q] for q in gate.qubits)
        for q in gate.qubits: last_layer[q] = layer
    depth = max(last_layer, default=0); total = len(circuit.gates); toffoli = counts["CCX"]
    return {"n_qubits": circuit.n_qubits, "n_input": circuit.n_input,
            "n_output": circuit.n_output, "n_ancilla": circuit.n_ancilla,
            "total_gates": total, "X": counts["X"], "CX": counts["CX"],
            "Toffoli_CCX": toffoli, "logical_depth_all_to_all": depth,
            "Toffoli_times_qubits": toffoli * circuit.n_qubits,
            "gate_count_times_qubits": total * circuit.n_qubits,
            "metric_note": "Profundidade lógica; ignora conectividade, decomposição física e correção de erros."}
