#!/usr/bin/env python3
"""Checa somente precondicoes racionais do lema da cauda; nao certifica RH."""
from fractions import Fraction as F

L = F(9, 20)
R = F(160)
pi_upper = F(22, 7)
pi_lower = F(333, 106)
assert pi_lower < pi_upper  # limites classicos estritos de pi

for n in (6, 8, 10, 20):
    assert R >= 1 and R > n*pi_upper/L, f"N={n} falhou"
assert R < 30*pi_lower/L, "N=30 foi aceito incorretamente"
# sinal(R_K)=(-1)^(K+1): K=8 -> negativo
assert (-1)**(8+1) == -1
print("PASS: N=6,8,10,20 admissiveis; N=30 nao em R=160; sinal K=8 negativo.")
print("AVISO: testes necessarios, nao certificacao integral nem RH.")
