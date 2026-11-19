#!/usr/bin/env python3
"""P246: rigorous scaling audit for the cofinal high-frequency cutoff.

For rational L > 0, compute
  A_L = 2 sum_{n < exp(2L)} Lambda(n)/sqrt(n)
with Arb balls, choose an exact rational T >= max(504, 4*pi*exp(A_L)),
and select N=ceil(L*T). The P239 Bessel-tail ratio condition is then checked
exactly with fractions for both parity sectors. This is a scaling audit only;
it does not certify positivity of the finite low-frequency block.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

from flint import arb, ctx

RATIO = F(3, 40)
MIN_T = F(504)


def aq(q: F | int) -> arb:
    q = F(q)
    return arb(q.numerator) / arb(q.denominator)


def dyadic(m: int, e: int) -> F:
    return F(m << e) if e >= 0 else F(m, 1 << -e)


def arb_upper_fraction(x: arb) -> F:
    """Return an exact rational upper bound for an Arb ball."""
    m, e = x.mid().man_exp()
    r, k = x.rad().upper().man_exp()
    return dyadic(int(m), int(e)) + dyadic(int(r), int(k))


def arb_wire(x: arb) -> dict[str, list[str] | str]:
    m, e = x.mid().man_exp()
    r, k = x.rad().upper().man_exp()
    return {
        "mid_dyadic": [str(m), str(e)],
        "radius_upper_dyadic": [str(r), str(k)],
        "display": str(x),
    }


def primes_up_to(nmax: int) -> list[int]:
    if nmax < 2:
        return []
    sieve = bytearray(b"\x01") * (nmax + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(nmax) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:nmax + 1:p] = b"\x00" * (((nmax - start) // p) + 1)
    return [i for i, is_prime in enumerate(sieve) if is_prime]


def cofinal_row(L: F, bits: int = 512) -> dict:
    if L <= 0:
        raise ValueError("L must be a positive rational")

    # Isolate floor(exp(2L)) by increasing Arb precision if necessary.
    exponent = 2 * L
    active_bits = bits
    while True:
        ctx.prec = active_bits
        X = aq(exponent).exp()
        floor_x = int(X.lower().floor().fmpz())
        floor_hi = int(X.upper().floor().fmpz())
        if floor_x == floor_hi and not X.contains_integer():
            break
        active_bits *= 2
        if active_bits > bits * 64:
            raise ArithmeticError("could not certify the strict prime-power cutoff")

    # For rational L>0, exp(2L) is not an integer. Thus n<exp(2L) is
    # exactly n<=floor(exp(2L)); sum Lambda(n) over its prime-power terms.
    A = arb(0)
    term_count = 0
    for p in primes_up_to(floor_x):
        logp = arb(p).log()
        n = p
        while n <= floor_x:
            A += 2 * logp / arb(n).sqrt()
            term_count += 1
            n *= p

    # Choose a rational strictly above 4*pi*exp(A_L), with relative slack
    # beyond the Arb upper endpoint so the inclusion remains decidable at scale.
    requested_T = 4 * arb.pi() * A.exp()
    requested_T_upper = arb_upper_fraction(requested_T)
    slack = max(F(1), requested_T_upper / (1 << (active_bits // 2)))
    T = max(MIN_T, requested_T_upper + slack)
    cutoff_gap = aq(T) - requested_T
    if not cutoff_gap.lower() > 0:
        raise ArithmeticError("rational cutoff failed to dominate 4*pi*exp(A_L)")

    x = L * T
    N = (x.numerator + x.denominator - 1) // x.denominator
    sectors = []
    for parity in (0, 1):
        n0 = 2 * N + parity
        ratio_sq = x**4 / ((2*n0 + 1) * (2*n0 + 3)**2 * (2*n0 + 5))
        if not ratio_sq < RATIO**2:
            raise ArithmeticError(f"Bessel ratio gate failed for parity {parity}")
        sectors.append({
            "parity": "even" if parity == 0 else "odd",
            "first_discarded_degree": str(n0),
            "ratio_squared_exact": str(ratio_sq),
            "ratio_below_3_over_40": True,
        })

    Tball = aq(T)
    log10_T = Tball.log() / arb(10).log()
    return {
        "L": f"{L.numerator}/{L.denominator}",
        "precision_bits_used": active_bits,
        "exp_2L_ball": arb_wire(X),
        "prime_power_cutoff_max_n": floor_x,
        "number_of_prime_powers_in_sum": term_count,
        "A_L_ball": arb_wire(A),
        "A_L_asymptotic_normalization": "A_L/(4*exp(L)) -> 1 by PNT and Abel summation",
        "T_rational_upper_cutoff": str(T),
        "T_cutoff_guarantee": "T >= 504 and T > 4*pi*exp(A_L), confirmed by Arb",
        "log10_T_ball": arb_wire(log10_T),
        "L_times_T_exact": str(x),
        "modes_per_parity_sufficient": str(N),
        "mode_cutoff_rule": "N=ceil(L*T); exact P239 ratio gate checked for both parities",
        "sector_ratio_checks": sectors,
        "scope": "scaling and tail-ratio gate only; no finite low-frequency block positivity or Weil-global positivity is certified",
    }


def parse_rational(text: str) -> F:
    try:
        return F(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError(f"invalid rational L: {text}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", action="append", type=parse_rational,
                        help="positive rational half-width; repeat to add rows")
    parser.add_argument("--bits", type=int, default=512)
    parser.add_argument("--out", type=Path, default=Path("numeria/P246_COFINAL_SCALING.json"))
    args = parser.parse_args()
    if args.bits < 128:
        parser.error("--bits must be at least 128")
    widths = args.L or [F(161, 200), F(1), F(3, 2), F(2), F(5, 2), F(3), F(4)]
    report = {
        "schema": "NUMERIA.ES8P246.COFINAL_SCALING.1",
        "arithmetic": "python-flint Arb for A_L and log10(T); exact Fraction for rational cutoffs and Bessel ratio checks",
        "prime_number_theorem_input": "psi(x)=sum_{n<=x}Lambda(n)~x",
        "derived_asymptotic": "A_L=2*sum_{n<exp(2L)}Lambda(n)/sqrt(n)~4*exp(L); log(T_L)=4*exp(L)+o(exp(L)) for T_L=4*pi*exp(A_L)",
        "mode_scale": "the checked sufficient choice N=ceil(L*T) has the same double-exponential scale as L*T; it says nothing about positivity of the low-frequency block",
        "rows": [cofinal_row(L, args.bits) for L in widths],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    for row in report["rows"]:
        print(f"L={row['L']}  cutoff_n<={row['prime_power_cutoff_max_n']}  "
              f"A_L={row['A_L_ball']['display']}  log10(T)={row['log10_T_ball']['display']}  "
              f"N/parity={row['modes_per_parity_sufficient']}")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
