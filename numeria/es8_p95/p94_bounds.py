#!/usr/bin/env python3
"""ES8-P95: exact rational bounds from ES8-P94 (outer Fourier quadrature).

A prior numerical error bound does not certify Arb moment evaluation.
"""
from fractions import Fraction as Q
from math import factorial

R = 600
SEGMENTS = 240
SEG_LEN = Q(5, 2)
H = Q(2, 5)
L = Q(5, 4)
K = 14
M = 1024

def panel_lengths_for_segment():
    return (Q(1, 5),)*12 + (Q(1, 10),)

def moment_taylor_error():
    return 2*H * (H*L)**(K+1) / factorial(K+1)

def panel_gl_error(length, n=12):
    if n < 1 or length <= 0 or length > Q(1, 5):
        raise ValueError("unproved GL panel")
    return (Q(factorial(n)**4, (2*n+1)*factorial(2*n)**2)
            * M * 4**(2*n) * length**(2*n+1))

def segment_gl_error(n=12):
    return sum((panel_gl_error(x, n) for x in panel_lengths_for_segment()), Q())

def full_gl_error(n=12):
    return SEGMENTS * segment_gl_error(n)

def proven_preconditions():
    assert sum(panel_lengths_for_segment(), Q()) == SEG_LEN
    assert len(panel_lengths_for_segment()) == 13
    assert Q(613*4, 3) < M
    assert Q(4, 5)/(1-Q(1,10)) < 1
    assert moment_taylor_error() < Q(2,10**17)
    assert full_gl_error(12) < Q(16,10**12)
    return {
        "status": "COTAS_RACIONAIS_ANALITICAS_VERIFICADAS",
        "scope": "outer GL 0..600 plus Taylor remainder for bump transform; not complete W12",
        "gl_total_error_upper_exact": str(full_gl_error(12)),
        "taylor_B_remainder_exact": str(moment_taylor_error()),
        "number_of_segments": SEGMENTS,
        "panels_per_fourier_segment": 13,
        "gauss_legendre_order": 12,
        "max_panel_length": "1/5",
        "majorant_M": M,
    }

if __name__ == "__main__":
    import json
    print(json.dumps(proven_preconditions(), indent=2))
