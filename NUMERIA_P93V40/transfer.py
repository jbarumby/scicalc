"""Directed Schur and strip transfer budgets, conditional on the complete heads."""

from fractions import Fraction
from functools import lru_cache
import json
from math import factorial
from pathlib import Path

from flint import acb, arb, ctx


PACKAGE = Path(__file__).resolve().parents[1] / "NUMERIA_P93V39"
T = 4500
BETA = 2
RHO = Fraction(33, 25)
WIDTH_CAP = Fraction(1387, 1000)
H = Fraction(1, 10**60)
TARGET = Fraction(9, 10) * H
THETAS = tuple(Fraction(j, 4) for j in range(5))
PRIMES = ((2, 2), (3, 3), (4, 2), (5, 5), (7, 7), (8, 2), (9, 3), (11, 11), (13, 13))


def ball(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator


def dyadic(x):
    m, e = x.man_exp()
    return [str(m), int(e)]


def count_mass():
    return sum((2 * arb(p).log() / arb(q).sqrt() for q, p in PRIMES), arb(0))


def low_symbol_cap(mass):
    """Use monotonicity of Re psi(1/4+iy) and the prime triangle bound."""
    quarter = arb(1) / 4
    pi_log = arb.pi().log()
    lo = quarter.digamma() - pi_log - BETA - mass
    hi = acb(quarter, arb(T) / 2).digamma().real - pi_log - BETA + mass
    cap = max((-lo).upper(), hi.upper())
    if not cap < 18:
        raise ArithmeticError("low spectral cap did not close")
    return cap


@lru_cache(maxsize=None)
def pole_audit(n):
    factorial_n = factorial(2 * n)
    return Fraction(48, factorial_n), Fraction(72, factorial_n**2)


def certify_transfer():
    with ctx.workprec(512):
        historical = json.loads((PACKAGE / "hybrid39.json").read_text())
        c, d, lam = (ball(historical[key]) for key in ("c", "d", "lambda"))
        f1 = ball(historical["f1"])
        width = Fraction(json.loads((PACKAGE / "moments_adaptive_binet_12.json").read_text())["width"])
        delta = 2 * ball(width) - arb(16).log()
        if not (delta > 0 and ball(width) < ball(WIDTH_CAP) and 2 * ball(width) < arb(17).log()):
            raise ArithmeticError("width does not cross precisely the threshold at 16")
        if not arb(2).log() / 4 < ball("1/5") or not arb.pi() > 3:
            raise ArithmeticError("strip bound constants did not close")
        if not (arb(2) * ball(WIDTH_CAP) * ball(WIDTH_CAP).exp() < 15 and
                arb(2) * ball(WIDTH_CAP).sinh() < 9):
            raise ArithmeticError("pole constants did not close")
        omega = 4_000_000
        mass = count_mass()
        # Discard the positive portion of the real Binet remainder. A full
        # absolute bound 1/(3*t) also gives this conservative lower bound.
        bfloor = (arb(omega) / (2 * arb.pi())).log() - 2 / arb(omega)**2 - 1 / (3 * arb(omega))
        if not bfloor - mass > BETA + 1:
            raise ArithmeticError("far spectral reserve did not close")
        loss = (omega * delta / 6).upper()
        transferred_margin = (ball(TARGET) - loss).lower()
        if not transferred_margin > 0:
            raise ArithmeticError("the conditional crossing margin is exhausted")

        rho = ball(RHO)
        a = ball(WIDTH_CAP)
        g_factor = 8 * a * T / 3 * (a * T * (rho - 1 / rho)).exp()
        g_factor /= (1 - 1 / rho)**2 * (1 - rho**-4)
        low_cap = low_symbol_cap(mass)

        def gate(n, theta, spectral_cap):
            gg = (g_factor * rho**(-4*n)).upper()  # Even parity is the worse bound.
            pole_cross, pole_tail = map(ball, pole_audit(n))
            mixing = spectral_cap + ball(theta) * d / (1 - lam)
            reserve = (BETA - (1-ball(theta))*c - ball(theta)*f1
                       - spectral_cap*gg - pole_tail - ball(TARGET))
            cross_square = 2 * mixing**2 * gg + 2 * pole_cross**2
            slack = (ball(H - TARGET)*reserve - cross_square).lower()
            return bool(slack > 0 and reserve > 0), gg, reserve.lower(), slack

        reports = []
        for theta in THETAS:
            row = {"theta": str(theta)}
            for cap in (50, 18):
                a_n, b_n = 3000, 3400
                if not gate(b_n, theta, cap)[0]:
                    raise ArithmeticError("the dimension bracket did not close")
                while b_n - a_n > 1:
                    mid = (a_n + b_n) // 2
                    if gate(mid, theta, cap)[0]:
                        b_n = mid
                    else:
                        a_n = mid
                _, g, reserve, slack = gate(b_n, theta, cap)
                row[f"symbol_cap_{cap}"] = {
                    "first_sufficient_dimension": b_n,
                    "g_even_upper_dyadic": dyadic(g),
                    "tail_reserve_lower_dyadic": dyadic(reserve),
                    "schur_slack_lower_dyadic": dyadic(slack),
                    "passes_at_3309": gate(3309, theta, cap)[0],
                }
            reports.append(row)
        return {
            "schema": "NUMERIA.P93V40.TRANSFER.1",
            "conditional_on_complete_heads_and_spectral_majorant": True,
            "width_exact": str(width),
            "rho_exact": str(RHO),
            "head_target_exact": str(H),
            "delta_lower_dyadic": dyadic(delta.lower()),
            "delta_upper_dyadic": dyadic(delta.upper()),
            "far_symbol_lower_dyadic": dyadic((bfloor-mass).lower()),
            "low_abs_symbol_upper_dyadic": dyadic(low_cap),
            "new_prime_loss_upper_dyadic": dyadic(loss),
            "conditional_crossing_margin_lower_dyadic": dyadic(transferred_margin),
            "theta_results": reports,
        }


if __name__ == "__main__":
    print(json.dumps(certify_transfer(), indent=2))
