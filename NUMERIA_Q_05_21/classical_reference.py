#!/usr/bin/env python3
"""Referência clássica independente: inverso modular via algoritmo euclidiano estendido."""
from __future__ import annotations

def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Retorna (g,x,y) com ax+by=g=gcd(a,b), usando iteração euclidiana."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t

def modular_inverse_euclid(x: int, p: int) -> int:
    """Calcula x^-1 mod p, levantando ValueError se o inverso não existir."""
    g, a, _ = extended_gcd(x, p)
    if g != 1:
        raise ValueError(f"Inverso inexistente: gcd({x},{p})={g}")
    return a % p
