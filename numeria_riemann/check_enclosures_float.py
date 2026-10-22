"""Checagem INDEPENDENTE e NÃO rigorosa dos enclosures exportados (numpy/scipy, double).

Recalcula a matriz projetada A = B^T (G + P) B por outro caminho, sem importar o
código com Arb: quadratura de Gauss-Legendre composta para a parte arquimediana
(digamma complexo do scipy), quadratura direta para os overlaps de primos e para
os vetores de restrição de polo. Depois confere se cada valor em ponto flutuante cai
dentro do enclosure exportado (com tolerância FLOAT_TOL, pois as entradas nulas por
estrutura saem ~1e-15 em double) e mede a distância ao ponto médio.

Não prova nada: serve para pegar erros de álgebra, sinal ou normalização, que um
enclosure rigoroso de uma conta errada não detectaria.

Uso: python -m numeria_riemann.check_enclosures_float RELATORIO.json [...]
"""
import json
import sys
from fractions import Fraction
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy.special import digamma

_X, _W = np.polynomial.legendre.leggauss(16)


def _gl(f, a, b, panels):
    edges = np.linspace(a, b, panels + 1)
    total = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = 0.5 * (hi - lo) * _X + 0.5 * (hi + lo)
        total += 0.5 * (hi - lo) * np.sum(_W * f(x))
    return total


def _sinc(x):
    safe = np.where(np.abs(x) < 1e-12, 1.0, x)
    return np.where(np.abs(x) < 1e-12, 1.0, np.sin(safe) / safe)


def _prime_power_log(n):
    m, d = n, 2
    while d * d <= m:
        if m % d == 0:
            while m % d == 0:
                m //= d
            return np.log(d) if m == 1 else None
        d += 1
    return np.log(n)


def projected_matrix(L, N, R, panels=6000):
    w = np.array([k * np.pi / (2 * L) for k in range(1, N + 1)])
    omega = lambda r: digamma(0.25 + 0.5j * r).real - np.log(np.pi)
    amp = lambda r, k: 2 * np.sqrt(L) * w[k] / (w[k] + r) * _sinc(L * (r - w[k]))
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            d = i - j
            if d % 2:
                continue
            phase = -1.0 if (abs(d) // 2) % 2 else 1.0
            G[i, j] = G[j, i] = _gl(lambda r: phase * amp(r, i) * amp(r, j) * omega(r) / np.pi, 0.0, R, panels)
    P = np.zeros((N, N))
    n = 2
    while np.log(n) <= 2 * L:
        lam = _prime_power_log(n)
        if lam is not None:
            a = np.log(n)
            for i in range(N):
                for j in range(N):
                    c = _gl(lambda x: np.sin(w[i] * (x + a)) * np.sin(w[j] * x), 0.0, 2 * L - a, 40) / L
                    P[i, j] += -2 * lam / np.sqrt(n) * c
        n += 1
    P = (P + P.T) / 2
    # restrições F_c(+-i/2) = 0, por quadratura direta (sem fórmula fechada)
    cons = np.array([[_gl(lambda x: np.sin(w[k] * x) / np.sqrt(L) * np.exp(s * (x - L) / 2), 0.0, 2 * L, 40)
                      for k in range(N)] for s in (+1.0, -1.0)])
    C2 = cons[:, :2]
    B = np.zeros((N, N - 2))
    for col in range(N - 2):
        B[:2, col] = np.linalg.solve(C2, -cons[:, 2 + col])
        B[2 + col, col] = 1.0
    assert np.allclose(cons @ B, 0.0, atol=1e-12)
    return B.T @ (G + P) @ B


def _float(text):
    return float(Fraction(Decimal(text)))


FLOAT_TOL = 1e-10  # erro esperado da quadratura em double (~1e-12) com folga; entradas nulas por estrutura saem ~1e-15


def check(report):
    p = report["parameters"]
    A = projected_matrix(float(Fraction(Decimal(str(p["L"])))), int(p["N"]), float(Fraction(Decimal(str(p["R"])))))
    rows = report["projected_matrix_enclosures"]
    n = len(rows)
    inside, worst_mid, worst_ratio = True, 0.0, 0.0
    for i in range(n):
        for j in range(n):
            lo, hi = _float(rows[i][j]["lower"]), _float(rows[i][j]["upper"])
            mid, rad = 0.5 * (lo + hi), 0.5 * (hi - lo)
            inside &= (lo - FLOAT_TOL) <= A[i, j] <= (hi + FLOAT_TOL)
            worst_mid = max(worst_mid, abs(A[i, j] - mid))
            worst_ratio = max(worst_ratio, abs(A[i, j] - mid) / max(rad, FLOAT_TOL))
    return {
        "float_tolerance": FLOAT_TOL,
        "all_float_values_inside_enclosures": bool(inside),
        "max_abs_distance_to_midpoint": worst_mid,
        "max_distance_over_radius": worst_ratio,
        "min_eigenvalue_float": float(np.linalg.eigvalsh((A + A.T) / 2).min()),
    }


def main(argv=None):
    ok = True
    for path in (argv or sys.argv[1:]):
        result = check(json.loads(Path(path).read_text(encoding="utf-8")))
        ok &= result["all_float_values_inside_enclosures"]
        print(path, json.dumps(result))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
