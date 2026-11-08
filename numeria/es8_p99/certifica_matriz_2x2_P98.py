#!/usr/bin/env python3
"""Enclosures Arb diretos de W11, W22, W12 e teste rigoroso do determinante.

Usa a implementação de integração P96 em W12_arb_direto_P96.py como motor
validado de quadraturas (python-flint==0.9.0). Calcula na normalização P92,
sem ler intervalos P89/P95, as entradas para os centros 1 e 2. A integração
usa a fórmula explícita P96 para cada par (n,m), a correlação compacta R(s),
o tratamento da origem e a cauda exata.

A decisão final usa apenas extremos racionais externos extraídos dos balls Arb:
  L_D = lower(W11)*lower(W22) - max(|lower(W12)|,|upper(W12)|)^2.
Se lower(W11)>0, lower(W22)>0 e L_D>0, a matriz simétrica 2x2 é positiva
 definida. Caso contrário, o resultado é INCONCLUSIVE (ou NOT_CERTIFIED), nunca
se inferindo positividade de autovalores float64.

Este script prova, no máximo, positividade no subespaço gerado pelas duas
funções de teste. Não prova positividade universal de Weil nem RH. A validade
matemática completa da fórmula explícita, integração aninhada e lemas de borda
continua sujeita a revisão independente, conforme advertência do motor P96.

Arquivos no mesmo diretório:
  W12_arb_direto_P96.py
Exemplo:
  python3 certifica_matriz_2x2_arb.py --prec 128 --tol 1e-6 \
      --eval-limit 12000 --depth-limit 20 --output matriz2x2.json
"""
from __future__ import annotations
import argparse
import json
import math
import sys
import time
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq
from W12_arb_direto_P98 import W12P92

H = Fraction(2, 5)
PRIME_POWERS = ((2, 2), (3, 3), (4, 2))  # (k, primo p), Lambda(k)=log(p)


def fmpq_fraction(q) -> Fraction:
    z = q.fmpq()
    return Fraction(int(z.p), int(z.q))


def endpoint_fraction(x: arb, upper: bool) -> Fraction:
    return fmpq_fraction(x.upper() if upper else x.lower())


def rational_decimal(q: Fraction, digits: int = 30, upper: bool = False) -> str:
    with localcontext() as dc:
        dc.prec = max(120, digits + 50)
        dc.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        scale = Decimal(1).scaleb(-digits)
        val = Decimal(q.numerator) / Decimal(q.denominator)
        return str(val.quantize(scale, rounding=dc.rounding))


def interval_record(z: acb) -> dict:
    if not z.is_finite():
        raise ArithmeticError(f"non-finite Arb result: {z}")
    if not z.imag.contains(0):
        raise ArithmeticError(f"real functional returned non-real enclosure: {z}")
    lo = endpoint_fraction(z.real, False)
    hi = endpoint_fraction(z.real, True)
    ilo = endpoint_fraction(z.imag, False)
    ihi = endpoint_fraction(z.imag, True)
    return {
        "lo": f"{lo.numerator}/{lo.denominator}",
        "hi": f"{hi.numerator}/{hi.denominator}",
        "im_lo": f"{ilo.numerator}/{ilo.denominator}",
        "im_hi": f"{ihi.numerator}/{ihi.denominator}",
        "lo_decimal_outward": rational_decimal(lo),
        "hi_decimal_outward": rational_decimal(hi, upper=True),
        "arb": str(z),
    }


def arb_bounds(rec: dict) -> tuple[Fraction, Fraction]:
    return Fraction(rec["lo"]), Fraction(rec["hi"])


def make_integrator(tol: str, eval_limit: int, depth_limit: int, eps_denominator: int = 2048) -> W12P92:
    return W12P92.make(tol, eval_limit, depth_limit, eps_denominator=eps_denominator)


def compute_entry(n: int, m: int, calc: W12P92) -> dict:
    """Compute one W_nm directly, in P92 normalization, returning Arb balls."""
    if not (n > 0 and m >= n):
        raise ValueError("centers must be positive integers with n<=m")

    h = calc.h
    delta_real = arb(m).log() - arb(n).log()
    delta = acb(delta_real)
    delta_abs = abs(delta_real)
    T = 2*h + delta_abs
    a = calc.moment_A()

    # Pole terms from the explicit formula.
    ratio = arb(m) / arb(n)
    pole = a*a*(ratio.sqrt() + (1/ratio).sqrt())

    # Prime powers k with Lambda(k)!=0. Correlation returns exact zero when
    # support proves that a term vanishes; including k=2,3,4 is sufficient for
    # h=2/5 and centers 1,2.
    prime_sum = acb(0)
    for k, p in PRIME_POWERS:
        logk = arb(k).log()
        coeff = arb(p).log()/arb(k).sqrt()
        prime_sum += coeff*(calc.correlation(acb(logk)+delta)
                            + calc.correlation(acb(logk)-delta))

    rdelta = calc.R_delta(delta)
    gamma_log = arb.const_euler() + (4*arb.pi()).log()

    # Stable origin split. The omitted correction is bounded by eps^2/h,
    # via the P96 rational lemma; the leading K contribution is integrated.
    eps = calc.eps
    small_main = rdelta*calc.small_kernel_integral()
    eps_fraction = Fraction(1, calc.eps_denominator)
    correction_bound = eps_fraction*eps_fraction/H
    correction_radius = arb(fmpq(correction_bound.numerator,
                                  correction_bound.denominator))
    small_correction = acb(arb(0, correction_radius), arb(0))

    # The first support contact is 2h-|delta|. For delta=0 this is T, so the
    # construction naturally has no duplicate zero/contact breakpoint.
    first_contact = 2*h-delta_abs
    if not bool(first_contact.lower() > eps.upper()):
        raise ArithmeticError("origin/contact panels overlap; refine eps or revise partition")
    points = calc.geometric_panels(eps, first_contact)

    # The t=|delta| breakpoint is optional for smooth R, but retaining it
    # separates the central correlation geometry when it lies after first_contact.
    if bool(delta_abs.lower() > points[-1].upper()) and bool(delta_abs.upper() < T.lower()):
        points.append(delta_abs)
    if bool(T.lower() > points[-1].upper()):
        points.append(T)
    if not all(bool(points[i].upper() < points[i+1].lower())
               for i in range(len(points)-1)):
        raise ArithmeticError("outer interval panels overlap or are unordered")

    I_mid = acb(0)
    for lo, hi in zip(points[:-1], points[1:]):
        I_mid += calc.integrate(calc.outer_integrand(delta, rdelta, lo, hi),
                                lo, hi, f"I_{n}{m} panel [{lo},{hi}]")

    # Exact compact-support tail: F_delta(t)=0 for t>=T.
    tail = rdelta*((T/2).tanh()).log()
    I_total = small_main + small_correction + I_mid + tail
    entry = pole - prime_sum - gamma_log*rdelta - I_total
    if not entry.is_finite():
        raise ArithmeticError(f"W_{n}{m} non-finite")

    return {
        "n": n, "m": m,
        "W": interval_record(entry),
        "components": {
            "moment_A": interval_record(a),
            "pole": interval_record(pole),
            "prime_sum": interval_record(prime_sum),
            "R_delta": interval_record(rdelta),
            "small_I_main": interval_record(small_main),
            "small_correction_radius_rational": f"{correction_bound.numerator}/{correction_bound.denominator}",
            "I_mid": interval_record(I_mid),
            "I_tail": interval_record(tail),
            "I_total": interval_record(I_total),
        },
        "delta": str(delta_real), "T": str(T),
        "integration_calls": calc.integral_count,
        "correlation_calls": calc.calls_R,
        "outer_panel_count": len(points)-1,
        "outer_analytic_denials": calc.outer_analytic_denials,
        "inner_support_trims": calc.inner_support_trims,
    }


def determinant_lower_bound(w11: dict, w22: dict, w12: dict) -> dict:
    l11,u11 = arb_bounds(w11); l22,u22 = arb_bounds(w22); l12,u12 = arb_bounds(w12)
    abs12 = max(abs(l12), abs(u12))
    d11_positive, d22_positive = l11 > 0, l22 > 0
    out = {
        "W11_lower_positive": d11_positive,
        "W22_lower_positive": d22_positive,
        "W12_abs_upper_rational": f"{abs12.numerator}/{abs12.denominator}",
        "derivation": "L_D = lower(W11)*lower(W22) - max(abs(lower(W12)),abs(upper(W12)))^2, exact rational arithmetic",
    }
    if not (d11_positive and d22_positive):
        out.update({
            "determinant_lower_rational": None,
            "determinant_lower_decimal_outward": None,
            "determinant_lower_strictly_positive": False,
            "reason": "diagonal positivity not certified; product of lower endpoints is not used as a determinant lower bound",
        })
        return out
    lower = l11*l22 - abs12*abs12
    out.update({
        "determinant_lower_rational": f"{lower.numerator}/{lower.denominator}",
        "determinant_lower_decimal_outward": rational_decimal(lower),
        "determinant_lower_strictly_positive": lower > 0,
    })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prec', type=int, default=128)
    ap.add_argument('--eps-den', type=int, default=2048)
    ap.add_argument('--tol', default='1e-6')
    ap.add_argument('--eval-limit', type=int, default=12000)
    ap.add_argument('--depth-limit', type=int, default=20)
    ap.add_argument('--output', type=Path, default=Path('W_2x2_ARB_P92.json'))
    args = ap.parse_args()
    if args.prec < 128 or args.eval_limit < 100 or args.depth_limit < 5:
        raise SystemExit('preconditions: prec>=128, eval_limit>=100, depth_limit>=5')
    ctx.prec = args.prec
    started = time.time()
    try:
        # Fresh state for each entry; no values from P89/P95 are read.
        entries = {}
        for key, n, m in (("W11",1,1),("W22",2,2),("W12",1,2)):
            calc = make_integrator(args.tol, args.eval_limit, args.depth_limit, args.eps_den)
            entries[key] = compute_entry(n,m,calc)
        det = determinant_lower_bound(entries['W11']['W'], entries['W22']['W'], entries['W12']['W'])
        positive = (det['W11_lower_positive'] and det['W22_lower_positive']
                    and det['determinant_lower_strictly_positive'])
        status = 'INTERVAL_POSITIVE_DEFINITE_CONDITIONAL' if positive else 'INCONCLUSIVE_DETERMINANT_NOT_POSITIVE'
        result = {
            "status": status,
            "version": "ES8-P99_P98_DIRECT_2X2_EXPERIMENT",
            "epsilon_denominator": args.eps_den,
            "precision_bits": args.prec,
            "tolerance": args.tol,
            "normalization": "b(x)=exp(-1/(1-(x/h)^2)); h=2/5",
            "centers": [1,2],
            "entries": entries,
            "determinant_test": det,
            "independence": "All entries computed directly from the P96 integral formula; no P89/P95 interval data loaded.",
            "scope": "At most certifies this finite 2x2 subspace. Does not establish Weil positivity for all test functions or RH.",
            "mathematical_review_warning": "The P96 engine itself flags the full Bombieri identity, nested integral parameterization and endpoint lemmas for independent mathematical review. Arb enclosures do not replace that review.",
            "elapsed_seconds": time.time()-started,
        }
    except Exception as exc:
        result = {
            "status": "INCONCLUSIVE",
            "reason": f"{type(exc).__name__}: {exc}",
            "elapsed_seconds": time.time()-started,
            "scope": "No positivity conclusion; an error or resource limit is not a certificate.",
        }
    text = json.dumps(result, indent=2, ensure_ascii=False)
    print(text)
    args.output.write_text(text+'\n', encoding='utf-8')
    print(f"Relatório salvo: {args.output}", file=sys.stderr)
    if result['status'] != 'INTERVAL_POSITIVE_DEFINITE_CONDITIONAL':
        raise SystemExit(2)

if __name__ == '__main__':
    main()