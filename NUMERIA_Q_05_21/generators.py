#!/usr/bin/env python3
"""Famílias independentes de síntese para inverso modular em campo primo."""
from __future__ import annotations
from numq_model import Circuit, Gate, append_mcx, append_pattern_controlled_x

def is_prime(n: int) -> bool:
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True

def dimensions(p: int) -> tuple[int, int]:
    if not is_prime(p): raise ValueError("O benchmark exige p primo")
    n = (p - 1).bit_length()
    return n, max(0, n - 2)

def inverse_table(p: int) -> tuple[int, list[int]]:
    n, _ = dimensions(p)
    table = [0] * (1 << n)
    for x in range(1, p): table[x] = pow(x, -1, p)
    return n, table

def _base_circuit(name: str, p: int):
    n, table = inverse_table(p)
    n_anc = max(0, n - 2)
    ins = tuple(range(n)); outs = tuple(range(n, 2*n))
    anc = tuple(range(2*n, 2*n+n_anc))
    return Circuit(name, n, n, n_anc, []), ins, outs, anc, table

def generate_minterm(p: int) -> Circuit:
    c, ins, outs, anc, table = _base_circuit("minterm_lookup", p)
    for x, y in enumerate(table):
        if not y: continue
        values = tuple((x >> i) & 1 for i in range(c.n_input))
        for bit, q_out in enumerate(outs):
            if (y >> bit) & 1:
                append_pattern_controlled_x(c.gates, ins, values, q_out, anc)
    return c

def anf_coefficients(values: list[int], n: int) -> list[int]:
    if len(values) != (1 << n): raise ValueError("Tabela-verdade incompatível com n")
    coeff = [int(v) & 1 for v in values]
    for bit in range(n):
        for mask in range(1 << n):
            if mask & (1 << bit): coeff[mask] ^= coeff[mask ^ (1 << bit)]
    return coeff

def generate_anf(p: int) -> Circuit:
    c, ins, outs, anc, table = _base_circuit("anf", p)
    for output_bit, q_out in enumerate(outs):
        truth = [(y >> output_bit) & 1 for y in table]
        coeff = anf_coefficients(truth, c.n_input)
        for mask, coefficient in enumerate(coeff):
            if not coefficient: continue
            controls = tuple(ins[i] for i in range(c.n_input) if (mask >> i) & 1)
            append_mcx(c.gates, controls, q_out, anc)
    return c
