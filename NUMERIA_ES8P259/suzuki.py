"""Certify Suzuki's scalar function on a finite interval in Arb."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from flint import arb, ctx


EPSILON = Fraction(1, 64)
PRECISION = 192
DEFAULT_LIMIT = 10**6
MAX_DEPTH = 24


def _ball(value):
    return arb(value.numerator) / value.denominator


def _dyadic(point):
    mantissa, exponent = point.man_exp()
    return [str(mantissa), int(exponent)]


def prime_powers(limit):
    """Return exactly the pairs (p^k,p) with 2 <= p^k <= limit."""
    if limit < 2:
        raise ValueError("the sieve limit must be at least 2")
    flags = bytearray(b"\1") * (limit + 1)
    flags[:2] = b"\0\0"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start::p] = b"\0" * ((limit - start) // p + 1)
    entries = []
    for p in range(2, limit + 1):
        if flags[p]:
            power = p
            while power <= limit:
                entries.append((power, p))
                power *= p
    return sorted(entries)


class SuzukiFormula:
    """Unconditional explicit formula in Theorem 1.1 of Suzuki."""

    def __init__(self):
        quarter = arb(1) / 4
        self.gamma_coefficient = (quarter.digamma() - arb.pi().log()) / 2
        self.lerch_at_zero = arb(2).zeta(quarter)
        self.tolerance = arb("1e-36")

    def initial_margin(self):
        """A strict rational lower bound for Psi(t)/t on (0, 1/64]."""
        if not self.gamma_coefficient > -arb(11) / 4:
            raise ArithmeticError("archimedean linear coefficient exceeds the bound")
        # (1 - exp(-lambda*t))/(lambda^2*t) >= 1/[lambda*(1+lambda*epsilon)].
        margin = sum((Fraction(256, (4 * n + 1) * (4 * n + 129)) for n in range(8)), Fraction(0))
        margin -= Fraction(11, 4)
        if not margin > 0 or not _ball(EPSILON) < arb(2).log():
            raise ArithmeticError("the initial no-prime interval failed")
        return margin

    def values(self, t):
        """Enclose the value and derivative of the prime-free part."""
        if not t > 0:
            raise ValueError("the Lerch series requires positive t")
        r = (-2 * t).exp()
        positive = (t / 2).exp()
        negative = (-t / 2).exp()
        geometric = 1 - r
        if not geometric > 0:
            raise ArithmeticError("the interval Lerch ratio does not contract")
        term = negative
        h1, h2 = arb(0), arb(0)
        for n in range(4097):
            reciprocal = arb(4) / (4 * n + 1)
            # Monotone reciprocal powers times geometric decay.
            error1 = (term * reciprocal / geometric).upper()
            if error1 < self.tolerance:
                error2 = (error1 * reciprocal).upper()
                h1 += arb(0, error1)
                h2 += arb(0, error2)
                break
            h1 += term * reciprocal
            h2 += term * reciprocal**2
            term *= r
        else:
            raise ArithmeticError("the Lerch tail did not reach the requested tolerance")
        value = 4 * (positive + negative - 2) + self.gamma_coefficient * t + (self.lerch_at_zero - h2) / 4
        slope = 2 * (positive - negative) + self.gamma_coefficient + h1 / 2
        return value, slope

    def tangent_lower(self, left, right, prime_slope, prime_offset, curvature_loss):
        """Use a lower tangent with a correction for possible negative curvature."""
        center, radius = (left + right) / 2, (right - left) / 2
        value, slope = self.values(center)
        value -= prime_slope * center - prime_offset
        slope -= prime_slope
        return (value - abs(slope) * radius - curvature_loss * radius**2 / 2).lower()


def certify(limit=DEFAULT_LIMIT):
    """Prove positivity on the continuous domain 0<t<=log(limit)."""
    powers = prime_powers(limit)
    with ctx.workprec(PRECISION):
        formula = SuzukiFormula()
        initial = formula.initial_margin()
        digest = sha256()
        min_location = None
        smallest = None
        count = depth_max = 0

        def cover(left, right, a, b, curvature, labels, code="", depth=0):
            nonlocal count, depth_max, smallest, min_location
            bound = formula.tangent_lower(left, right, a, b, curvature)
            if bound > 0:
                record = [*labels, code, _dyadic(bound)]
                digest.update((json.dumps(record, separators=(",", ":")) + "\n").encode())
                if smallest is None or bound < smallest:
                    smallest, min_location = bound, record[:3]
                count += 1
                depth_max = max(depth, depth_max)
                return
            if depth >= MAX_DEPTH:
                raise ArithmeticError(f"no positive bound for cell {labels}, subdivision {code}")
            center = (left + right) / 2
            cover(left, center, a, b, curvature, labels, code + "0", depth + 1)
            cover(center, right, a, b, curvature, labels, code + "1", depth + 1)

        log2 = arb(2).log()
        cover(_ball(EPSILON), log2, arb(0), arb(0), 33, ["1/64", "log(2)"])
        a, b = arb(0), arb(0)
        left = log2
        for index, (q, p) in enumerate(powers):
            if q == limit:
                break
            weight = arb(p).log() / arb(q).sqrt()
            a += weight
            b += weight * left
            right_q = powers[index + 1][0] if index + 1 < len(powers) else limit
            right = arb(right_q).log()
            if not left < right:
                raise ArithmeticError("prime power breakpoints are not separated")
            cover(left, right, a, b, 0, [f"log({q})", f"log({right_q})"])
            left = right
        return {
            "schema": "NUMERIA.ES8P259.SUZUKI_FINITE.1",
            "formula": "10.1112/jlms.12785, equations (1.1)-(1.3), Theorem 1.7",
            "precision_bits": PRECISION,
            "certified_initial_slope_lower_exact": str(initial),
            "initial_interval_end_exact": str(EPSILON),
            "positive_for_all_t_in": f"(0,log({limit})]",
            "prime_power_count_through_limit": len(powers),
            "certified_cells_after_initial_interval": count,
            "maximum_dyadic_subdivision_depth": depth_max,
            "minimum_verified_lower_dyadic_after_initial_interval": _dyadic(smallest),
            "minimum_bound_cell": min_location,
            "ordered_leaf_certificate_sha256": digest.hexdigest(),
            "unbounded_range_established": False,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = certify(args.limit)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
