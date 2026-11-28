"""Directed Schur and strip transfer budgets, conditional on the complete heads."""

from fractions import Fraction
from functools import lru_cache
import argparse
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


def certify_transfer(head_margin=H, *, separated_strips=False, chebyshev_tail=False,
                     retention=Fraction(9, 10)):
    head_margin = Fraction(head_margin)
    retention = Fraction(retention)
    if head_margin <= 0 or not 0 < retention < 1:
        raise ValueError("require a positive margin and a retention fraction in (0, 1)")
    target = retention * head_margin
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
        omega = 1_600_000 if separated_strips else 4_000_000
        mass = count_mass()
        # Discard the positive portion of the real Binet remainder. A full
        # absolute bound 1/(3*t) also gives this conservative lower bound.
        bfloor = (arb(omega) / (2 * arb.pi())).log() - 2 / arb(omega)**2 - 1 / (3 * arb(omega))
        if separated_strips:
            kappa = arb(2).log() / 4
            eta = ball("8/15")
            if not bfloor - mass - BETA > kappa * (1 + 1 / eta):
                raise ArithmeticError("far reserve does not absorb the high part on the strips")
            separation = arb(16).log() - delta
            if not separation > 0:
                raise ArithmeticError("positive separation of the strips is required")
            strip_bound = kappa * (1 + eta) * delta / arb.pi() * (omega + 1 / separation)
            trace_weight = (2 * kappa * (1 + eta) * delta).upper()
            trace_remainder = (4 * kappa * (1 + eta) * delta**3 * omega**3 /
                               (9 * arb.pi())).upper()
            trace_report = {
                "edge_evaluation_weight_upper_dyadic": dyadic(trace_weight),
                "edge_evaluation_scalar_remainder_upper_dyadic": dyadic(trace_remainder),
            }
        else:
            if not bfloor - mass > BETA + 1:
                raise ArithmeticError("far spectral reserve did not close")
            strip_bound = omega * delta / 6
        loss = strip_bound.upper()
        transferred_margin = (ball(target) - loss).lower()
        if not transferred_margin > 0:
            raise ArithmeticError("the conditional crossing margin is exhausted")

        rho = ball(RHO)
        d_rho = rho - 1 / rho
        if chebyshev_tail:
            # Cauchy's bound averaged on the ellipse, with parity and the
            # Chebyshev L² norm. Integrate separately below and above T/2.
            z = ball(width) * T * d_rho
            half = (z / 2).exp()
            g_factor = 2 * (half - 1 + arb.pi() * (z.exp() - half) / (2 * z))
            g_factor /= d_rho * (1 - rho**-4)
        else:
            a = ball(WIDTH_CAP)
            g_factor = 8 * a * T / 3 * (a * T * d_rho).exp()
            g_factor /= (1 - 1 / rho)**2 * (1 - rho**-4)
        low_cap = low_symbol_cap(mass)

        def gate(n, theta, spectral_cap):
            gg = (g_factor * rho**(-4*n)).upper()  # Even parity is the worse bound.
            pole_cross, pole_tail = map(ball, pole_audit(n))
            mixing = spectral_cap + ball(theta) * d / (1 - lam)
            reserve = (BETA - (1-ball(theta))*c - ball(theta)*f1
                       - spectral_cap*gg - pole_tail - ball(target))
            if chebyshev_tail:
                cross_square = (mixing * gg.sqrt() + pole_cross)**2
            else:
                cross_square = 2 * mixing**2 * gg + 2 * pole_cross**2
            slack = (ball(head_margin - target)*reserve - cross_square).lower()
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
            "head_target_exact": str(head_margin),
            "delta_lower_dyadic": dyadic(delta.lower()),
            "delta_upper_dyadic": dyadic(delta.upper()),
            "far_symbol_lower_dyadic": dyadic((bfloor-mass).lower()),
            "low_abs_symbol_upper_dyadic": dyadic(low_cap),
            "new_prime_loss_upper_dyadic": dyadic(loss),
            "conditional_crossing_margin_lower_dyadic": dyadic(transferred_margin),
            "theta_results": reports,
            **({"strip_method": "separate-diagonal-and-off-diagonal",
                "strip_cutoff": omega, "strip_young_exact": "8/15",
                "strip_bound_formula_lower_dyadic": dyadic(strip_bound.lower()),
                "edge_evaluation_alternative": trace_report}
               if separated_strips else {}),
            **({"chebyshev_tail": True} if chebyshev_tail else {}),
            **({"retention_exact": str(retention)} if retention != Fraction(9, 10) else {}),
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--margin", default=str(H))
    parser.add_argument("--separated-strips", action="store_true")
    parser.add_argument("--chebyshev-tail", action="store_true")
    parser.add_argument("--retain", type=Fraction, default=Fraction(9, 10))
    args = parser.parse_args()
    print(json.dumps(certify_transfer(args.margin, separated_strips=args.separated_strips,
                                      chebyshev_tail=args.chebyshev_tail, retention=args.retain),
                     indent=2))
