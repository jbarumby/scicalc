#!/usr/bin/env python3
"""Gauss-Legendre composta com resto clássico validado por Arb.

Este módulo implementa o mecanismo rigoroso da quadratura. O usuário deve
fornecer duas funções intervalares:
  f_interval(X): bola Arb que contém f(x) para TODO x em X;
  derivative_abs_bound(I): bola Arb não negativa cujo limite superior é
      >= sup_{x em I} |f^(2n)(x)|.

A rotina calcula raízes/pesos de Gauss-Legendre em bolas Arb (isolamento dos
zeros certificado por intervalos racionais com troca de sinal) e soma o resto
clássico em cada painel:
 |E_n| <= (n!)^4 / ((2n+1)((2n)!)^3) * length^(2n+1) * sup|f^(2n)|.

IMPORTANTE: sem uma implementação comprovadamente válida de derivative_abs_bound
para o integrando concreto, o resultado é apenas um enclosure condicional, não
uma certificação de W_12. O exemplo no fim valida a infraestrutura em f(x)=x^2.

Instalação: python -m pip install python-flint
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Sequence

from flint import arb, fmpq, ctx


def qarb(q: Fraction) -> arb:
    return arb(fmpq(q.numerator, q.denominator))


def qball(a: Fraction, b: Fraction) -> arb:
    if a > b:
        raise ValueError("intervalo racional invertido")
    midpoint = (a+b)/2
    radius = (b-a)/2
    return arb(fmpq(midpoint.numerator, midpoint.denominator),
               fmpq(radius.numerator, radius.denominator))


def upper_fraction(x: arb) -> Fraction:
    q = x.upper().upper().fmpq()
    return Fraction(int(q.p), int(q.q))


def lower_fraction(x: arb) -> Fraction:
    q = x.lower().lower().fmpq()
    return Fraction(int(q.p), int(q.q))


def fraction_text(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def legendre_coefficients(n: int) -> list[Fraction]:
    """Coeficientes racionais de P_n(x), em ordem crescente de potência."""
    if n < 0:
        raise ValueError("n precisa ser >= 0")
    p0 = [Fraction(1)]
    if n == 0:
        return p0
    p1 = [Fraction(0), Fraction(1)]
    if n == 1:
        return p1
    for k in range(1, n):
        # P_{k+1}=((2k+1)x P_k - k P_{k-1})/(k+1)
        xp = [Fraction(0)] + p1
        size = max(len(xp), len(p0))
        xp += [Fraction(0)] * (size-len(xp))
        p0pad = p0 + [Fraction(0)] * (size-len(p0))
        p2 = [((2*k+1)*xp[i]-k*p0pad[i])/Fraction(k+1) for i in range(size)]
        p0, p1 = p1, p2
    return p1


def eval_poly_q(coeff: Sequence[Fraction], x: Fraction) -> Fraction:
    y = Fraction(0)
    for a in reversed(coeff):
        y = y*x+a
    return y


def eval_poly_arb(coeff: Sequence[Fraction], x: arb) -> arb:
    y = arb(0)
    for a in reversed(coeff):
        y = y*x+qarb(a)
    return y


def derivative_coefficients(coeff: Sequence[Fraction]) -> list[Fraction]:
    return [Fraction(k)*coeff[k] for k in range(1, len(coeff))]


def isolate_legendre_roots(n: int) -> list[arb]:
    """Isola todos los n ceros simples de P_n con intervalos racionales.

    NumPy solo propone semillas aproximadas; la validez se basa en verificar,
    con aritmética racional exacta, cambio de signo en extremos disjuntos.
    n intervalos disjuntos con cambio de signo prueban que contienen los n
    ceros de un polinomio de grado n, uno por intervalo.
    """
    if n < 1:
        raise ValueError("la regla de Gauss necesita n>=1")
    try:
        import numpy as np
        seeds = np.polynomial.legendre.leggauss(n)[0]
    except Exception as exc:
        raise RuntimeError("NumPy se usa solo para semillas de raíces; instale numpy") from exc
    coeff = legendre_coefficients(n)
    roots: list[arb] = []
    brackets: list[tuple[Fraction, Fraction]] = []
    for seed in sorted(float(v) for v in seeds):
        center = Fraction(format(seed, ".17g"))
        width = Fraction(1, 10**14)
        found = None
        for _ in range(100):
            a, b = center-width, center+width
            fa, fb = eval_poly_q(coeff, a), eval_poly_q(coeff, b)
            if fa == 0:
                # Raíz racional exacta; crear bracket degenerado.
                found = (a, a)
                break
            if fb == 0:
                found = (b, b)
                break
            if (fa < 0 < fb) or (fb < 0 < fa):
                found = (a, b)
                break
            width *= 2
        if found is None:
            raise RuntimeError(f"não foi possível isolar a raiz {seed}")
        # Refina o bracket com bisseção racional exata. As sementes float não
        # entram na prova; só foram usadas para encontrar o primeiro bracket.
        a, b = found
        if a != b:
            fa = eval_poly_q(coeff, a)
            while b-a > Fraction(1, 10**40):
                mid = (a+b)/2
                fm = eval_poly_q(coeff, mid)
                if fm == 0:
                    a = b = mid
                    break
                if (fa < 0 < fm) or (fm < 0 < fa):
                    b = mid
                else:
                    a = mid
                    fa = fm
            found = (a,b)
        brackets.append(found)
    brackets.sort()
    for i, (a,b) in enumerate(brackets):
        if i and brackets[i-1][1] >= a:
            raise RuntimeError("intervalos das raízes não são disjuntos")
        if a == b:
            roots.append(qarb(a))
        else:
            # Prova de troca de sinal já feita exatamente nos racionais.
            roots.append(qball(a,b))
    if len(roots) != n:
        raise RuntimeError("número inesperado de raízes isoladas")
    return roots


def gauss_legendre_rule_arb(n: int) -> tuple[list[arb], list[arb]]:
    """Enclosures Arb dos nós e pesos da regra de Gauss-Legendre em [-1,1]."""
    p = legendre_coefficients(n)
    dp = derivative_coefficients(p)
    xs = isolate_legendre_roots(n)
    ws = []
    for x in xs:
        deriv = eval_poly_arb(dp, x)
        den = (1-x*x)*deriv*deriv
        if not (den.lower() > 0):
            raise RuntimeError("peso GL não ficou positivo; aumente ctx.prec")
        ws.append(2/den)
    return xs, ws


def gl_remainder_constant(n: int) -> Fraction:
    """C_n=(n!)^4/((2n+1)((2n)!)^3), exato racional."""
    return Fraction(math.factorial(n)**4,
                    (2*n+1)*math.factorial(2*n)**3)


@dataclass
class PanelResult:
    a: Fraction
    b: Fraction
    quadrature: arb
    remainder_radius: Fraction
    derivative_bound_upper: Fraction
    result: arb


def composite_gauss_legendre(
    f_interval: Callable[[arb], arb],
    derivative_abs_bound: Callable[[Fraction, Fraction], arb],
    endpoints: Sequence[Fraction],
    order: int,
) -> tuple[arb, list[PanelResult]]:
    """Integra em painéis racionais; retorno inclui soma e restos rigorosos.

    derivative_abs_bound(a,b) deve fornecer uma bola não negativa com limite
    superior M >= sup_[a,b]|f^(2*order)|. A responsabilidade por essa prova é
    de quem fornece a função; não é inferida automaticamente pelo Arb.
    """
    if order < 1:
        raise ValueError("order deve ser >=1")
    if len(endpoints) < 2 or any(endpoints[k] >= endpoints[k+1]
                                for k in range(len(endpoints)-1)):
        raise ValueError("endpoints devem ser estritamente crescentes")
    xs, ws = gauss_legendre_rule_arb(order)
    C = qarb(gl_remainder_constant(order))
    panel_results: list[PanelResult] = []
    total = arb(0)
    for a,b in zip(endpoints[:-1], endpoints[1:]):
        length = b-a
        mid = (a+b)/2
        half = length/2
        mid_ball, half_ball = qarb(mid), qarb(half)
        qsum = arb(0)
        for x,w in zip(xs,ws):
            X = mid_ball + half_ball*x
            fx = f_interval(X)
            qsum += half_ball*w*fx
        Mball = derivative_abs_bound(a,b)
        M = upper_fraction(Mball.upper())
        if M < 0:
            raise ValueError("cota de derivada negativa")
        E = gl_remainder_constant(order) * length**(2*order+1) * M
        # A adição de [-E,E] inclui rigorosamente o erro do teorema clássico.
        enclosure = qsum + arb(0, fmpq(E.numerator, E.denominator))
        panel_results.append(PanelResult(a,b,qsum,E,M,enclosure))
        total += enclosure
    return total, panel_results


def arb_interval_json(x: arb) -> dict:
    return {"lo": fraction_text(lower_fraction(x)),
            "hi": fraction_text(upper_fraction(x)),
            "arb": str(x)}


def demo() -> None:
    """Demonstração de infraestrutura com integral x^2 em [-1,1]."""
    # Para ordem 2, f^(4)=0 exatamente, logo o resto clássico é zero.
    def f_iv(X: arb) -> arb:
        return X*X
    def d_bound(a: Fraction, b: Fraction) -> arb:
        return arb(0)
    ans, panels = composite_gauss_legendre(
        f_iv, d_bound, [Fraction(-1), Fraction(1)], order=2)
    truth = arb(2)/3
    print(json.dumps({
        "status": "DEMO_ONLY",
        "integral_x2_minus1_1": arb_interval_json(ans),
        "contains_exact_2_over_3": bool(ans.contains(truth)),
        "panels": [{"a":fraction_text(p.a),"b":fraction_text(p.b),
                    "quadrature":str(p.quadrature),
                    "remainder_radius":fraction_text(p.remainder_radius)}
                   for p in panels],
        "warning": "Demonstra apenas a regra genérica. Não implementa a cota de derivadas do integrando W12."
    }, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--demo", action="store_true", help="testa integral de x^2")
    args = ap.parse_args()
    ctx.prec = args.prec
    if args.demo:
        demo()
    else:
        ap.print_help()
        print("\nAPI Python: composite_gauss_legendre(f_interval, derivative_abs_bound, endpoints, order)")


if __name__ == "__main__":
    main()