#!/usr/bin/env python3
"""Modelo mínimo de circuitos reversíveis clássicos sobre base computacional."""
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class Gate:
    name: str
    qubits: tuple[int, ...]
    def label(self) -> str:
        return f"{self.name}({','.join(map(str, self.qubits))})"

@dataclass
class Circuit:
    name: str
    n_input: int
    n_output: int
    n_ancilla: int
    gates: list[Gate]
    @property
    def n_qubits(self) -> int:
        return self.n_input + self.n_output + self.n_ancilla

def validate_gate(g: Gate, n_qubits: int) -> None:
    expected = {"X": 1, "CX": 2, "CCX": 3}
    if g.name not in expected or len(g.qubits) != expected[g.name]:
        raise ValueError(f"Porta inválida: {g}")
    if len(set(g.qubits)) != len(g.qubits):
        raise ValueError(f"Um gate não pode repetir qubit: {g}")
    if any(q < 0 or q >= n_qubits for q in g.qubits):
        raise ValueError(f"Qubit fora do circuito n={n_qubits}: {g}")

def apply_gate(state: int, g: Gate) -> int:
    """Aplica X, CX ou CCX a um vetor de bits; q0 é o bit menos significativo."""
    if g.name == "X":
        return state ^ (1 << g.qubits[0])
    if g.name == "CX":
        c, t = g.qubits
        return state ^ ((((state >> c) & 1)) << t)
    if g.name == "CCX":
        c1, c2, t = g.qubits
        return state ^ (((((state >> c1) & 1) & ((state >> c2) & 1))) << t)
    raise ValueError(f"Porta desconhecida: {g.name}")

def run_circuit(state: int, circuit: Circuit) -> int:
    for g in circuit.gates:
        validate_gate(g, circuit.n_qubits)
        state = apply_gate(state, g)
    return state

def append_mcx(gates: list[Gate], controls: tuple[int, ...], target: int,
               scratch: tuple[int, ...]) -> None:
    """Decompõe X com k controles em X/CX/CCX, limpando ancillas."""
    k = len(controls)
    if target in controls or len(set(controls)) != k:
        raise ValueError("Controles e alvo devem ser distintos")
    if k == 0:
        gates.append(Gate("X", (target,)))
    elif k == 1:
        gates.append(Gate("CX", (controls[0], target)))
    elif k == 2:
        gates.append(Gate("CCX", (controls[0], controls[1], target)))
    else:
        need = k - 2
        if len(scratch) < need or any(q in controls or q == target for q in scratch[:need]):
            raise ValueError(f"MCX com {k} controles precisa de {need} ancillas distintas")
        a = scratch[:need]
        gates.append(Gate("CCX", (controls[0], controls[1], a[0])))
        for j in range(2, k - 1):
            gates.append(Gate("CCX", (a[j - 2], controls[j], a[j - 1])))
        gates.append(Gate("CCX", (a[-1], controls[-1], target)))
        for j in range(k - 2, 1, -1):
            gates.append(Gate("CCX", (a[j - 2], controls[j], a[j - 1])))
        gates.append(Gate("CCX", (controls[0], controls[1], a[0])))

def append_pattern_controlled_x(gates: list[Gate], controls: tuple[int, ...],
                                values: tuple[int, ...], target: int,
                                scratch: tuple[int, ...]) -> None:
    """Aplica X no alvo quando controles==values (valores 0/1)."""
    if len(controls) != len(values):
        raise ValueError("Tamanhos diferentes de controles e padrão")
    negative = [q for q, v in zip(controls, values) if v == 0]
    for q in negative:
        gates.append(Gate("X", (q,)))
    append_mcx(gates, controls, target, scratch)
    for q in reversed(negative):
        gates.append(Gate("X", (q,)))
