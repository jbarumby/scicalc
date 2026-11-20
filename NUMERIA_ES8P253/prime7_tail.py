"""Certify the signed frequency tail with prime 7 and the P250 beta."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from flint import acb, arb, ctx

from NUMERIA_ES8P249.interval_audit import _arb_prime_terms, _dyadic
from NUMERIA_ES8P249.prime_operator import P236_PRIME_POWERS
from NUMERIA_ES8P250.adaptive_signed_tail import BETA_CUTOFF


POWERS = (*P236_PRIME_POWERS, (7, 7))
CUTOFF = 1551
END = 2200
WITNESS = Fraction(15501, 10)
PRECISION_BITS = 192
MAX_DEPTH = 20


def _rational(value):
    return arb(value.numerator) / value.denominator


def _binet_lower(frequency, pi):
    return (frequency / (2 * pi)).log() - 1 / frequency


def _positive_cover(terms, beta, cutoff, end, pi):
    curvature = sum(amplitude * shift**2 for amplitude, shift in terms)
    mass = sum(amplitude for amplitude, _ in terms)
    mass_lower = _binet_lower(arb(end), pi) - beta - mass
    if not mass_lower > 0:
        raise ArithmeticError("the monotone far-tail bound could not be certified")

    leaves = []
    attempts = 0
    deepest = 0

    def cover(left, right, depth):
        nonlocal attempts, deepest
        attempts += 1
        center = (left + right) / 2
        c, radius = _rational(center), _rational((right - left) / 2)
        oscillation, slope = arb(0), arb(0)
        for amplitude, shift in terms:
            sine, cosine = (c * shift).sin_cos()
            oscillation += amplitude * cosine
            slope -= amplitude * shift * sine
        lower = (
            _binet_lower(_rational(left), pi) - beta - oscillation
            - radius * abs(slope) - radius**2 * curvature / 2
        ).lower()
        if lower > 0:
            leaves.append((left, right, lower))
            deepest = max(deepest, depth)
            return
        if depth >= MAX_DEPTH:
            raise ArithmeticError(f"sign undecided on [{left}, {right}]")
        cover(left, center, depth + 1)
        cover(center, right, depth + 1)

    for start in range(cutoff, end):
        cover(Fraction(start), Fraction(start + 1), 0)
    return {
        "far_tail_start": end,
        "far_tail_uniform_strict_lower_dyadic": _dyadic(mass_lower.lower()),
        "adaptive_range": [cutoff, end],
        "adaptive_leaf_count": len(leaves),
        "adaptive_attempts": attempts,
        "adaptive_max_depth": deepest,
        "adaptive_uniform_strict_lower_dyadic": _dyadic(min(lower for _, _, lower in leaves)),
        "adaptive_cover": [[str(left), str(right), _dyadic(lower)] for left, right, lower in leaves],
    }


def certify_tail():
    """Prove the new cutoff and enclose a negative witness below it."""
    with ctx.workprec(PRECISION_BITS):
        pi = arb.pi()
        terms = _arb_prime_terms(POWERS)
        old_mass = sum(amplitude for amplitude, _ in terms[:-1])
        beta = _binet_lower(arb(BETA_CUTOFF), pi) - old_mass
        if not beta > 0:
            raise ArithmeticError("positive fixed beta could not be certified")
        point = _rational(WITNESS)
        oscillation = sum(amplitude * (point * shift).cos() for amplitude, shift in terms)
        witness = acb(arb(1) / 4, point / 2).digamma().real - pi.log() - beta - oscillation
        if not witness < 0:
            raise ArithmeticError("the exact-frequency negative witness could not be certified")

        return {
            "schema": "NUMERIA.ES8P253.PRIME7_TAIL.1",
            "scope": "sign of the Fourier symbol residual only",
            "precision_bits": PRECISION_BITS,
            "prime_powers": [q for q, _ in POWERS],
            "window_condition": "7 < exp(2L) < 8",
            "old_beta_reference_cutoff": BETA_CUTOFF,
            "old_beta_dyadic_endpoints": [_dyadic(beta.lower()), _dyadic(beta.upper())],
            "binet_bound": "Re digamma(1/4+it/2)-log(pi) >= log(t/(2pi))-1/t for t>0",
            "cutoff_minimal_integer_for_fixed_beta": CUTOFF,
            "certified_range_positive_frequency": [CUTOFF, "infinity"],
            "rational_negative_witness_frequency": str(WITNESS),
            "negative_witness_dyadic_endpoints": [_dyadic(witness.lower()), _dyadic(witness.upper())],
            **_positive_cover(terms, beta, CUTOFF, END, pi),
        }


def certify_four_prime_floor():
    """Certify the larger four-prime tail floor beta = 1/4 for use at cut 504."""
    with ctx.workprec(PRECISION_BITS):
        return {
            "schema": "NUMERIA.ES8P253.FOUR_PRIME_RAISED_FLOOR.1",
            "scope": "four-prime Fourier multiplier on the frequency tail only",
            "precision_bits": PRECISION_BITS,
            "prime_powers": [q for q, _ in P236_PRIME_POWERS],
            "proposed_beta_exact": "1/4",
            "certified_range_positive_frequency": [BETA_CUTOFF, "infinity"],
            **_positive_cover(_arb_prime_terms(P236_PRIME_POWERS), arb(1) / 4, BETA_CUTOFF, 650, arb.pi()),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--four-prime-floor", action="store_true")
    args = parser.parse_args()
    certificate = certify_four_prime_floor() if args.four_prime_floor else certify_tail()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(certificate, separators=(",", ":")) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
