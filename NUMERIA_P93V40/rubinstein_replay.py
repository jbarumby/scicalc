"""Exact checks of elementary inequalities from Rubinstein's zeta draft, v1."""

import argparse
from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path

from flint import arb, ctx

from .interval_ldl import dyadic


SOURCE = {
    "url": "https://math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper.pdf",
    "version": 1,
    "date": "2026-11-21",
    "pdf_pages": 164,
    "sha256": "667d88b3c375adf0466678094f452ca966aef0f4d368f6340d2bd5328b7fb706",
}


def check_top_radius():
    """Check (III.92) and (III.93) with a finite rational grid."""
    cases = 0
    for prime in (3, 5, 7):
        for cutoff in (1, 3, 5, 7, 9):
            exponents = tuple(range(1, cutoff, 2))
            lengths = (0, Fraction(1, 2), 1, prime, prime**3, prime**6)
            for length in map(Fraction, lengths):
                thresholds = {j: length / prime**j for j in exponents}
                for factors in product((0, Fraction(1, 2), 1, 2), repeat=len(exponents)):
                    capacities = {j: factor * thresholds[j]
                                  for j, factor in zip(exponents, factors, strict=True)}
                    first = next((j for j in exponents if capacities[j] >= thresholds[j]), cutoff)
                    radius = Fraction(1, prime**first)
                    cost = sum(capacities[j] * max(0, 1 - prime**j * radius) for j in exponents)
                    bound = sum(min(capacities[j], thresholds[j]) for j in exponents)
                    gamma = Fraction(prime, prime**2 - 1)
                    gain = sum(max(0, thresholds[j] - capacities[j]) for j in exponents)
                    gain += length * gamma / prime**(cutoff + 1)
                    if (radius < Fraction(1, prime**cutoff)
                            or length * radius + cost > bound + length / prime**cutoff
                            or length * (gamma - radius) - cost < gain):
                        raise ArithmeticError("top-radius inequality failed on the rational grid")
                    cases += 1
    return cases


def orbit_two_max_flow(prime, d0, d1):
    """Compute the LP maximum by enumerating intersections of its boundaries."""
    boundaries = ((1, 0, 0), (0, 1, 0), (prime, 1, d0), (1, prime, d1))
    values = []
    for (a, b, c), (d, e, f) in combinations(boundaries, 2):
        determinant = a * e - b * d
        if determinant == 0:
            continue
        x = Fraction(c * e - b * f, determinant)
        y = Fraction(a * f - c * d, determinant)
        if x >= 0 and y >= 0 and prime * x + y <= d0 and x + prime * y <= d1:
            values.append(x + y)
    return max(values)


def check_orbit_two():
    """Compare exact vertex enumeration with the closed form (III.97)."""
    cases = 0
    for prime in (3, 5, 7):
        for d0, d1 in product(range(33), repeat=2):
            if orbit_two_max_flow(prime, d0, d1) != min(d0, d1, Fraction(d0 + d1, prime + 1)):
                raise ArithmeticError("closed form disagrees with the rational LP optimum")
            cases += 1
    return cases


def check_prime_power_sums():
    """Recompute the four values of G(x) in §I.1 with an integer sieve and Arb."""
    targets = (10_000, 100_000, 1_000_000, 2_000_000)
    reported_millis = (1749, 3671, 11040, -12268)
    limit = targets[-1]
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = bytes(2)
    for prime in range(2, isqrt(limit) + 1):
        if flags[prime]:
            start = prime * prime
            flags[start::prime] = bytes((limit - start) // prime + 1)

    bits = 128
    with ctx.workprec(bits):
        sums = [arb(0) for _ in targets]
        reciprocal_sums = [arb(0) for _ in targets]
        counts = [0 for _ in targets]
        for prime, is_prime in enumerate(flags):
            if not is_prime:
                continue
            log_prime = arb(prime).log()
            power = prime
            while power <= limit:
                index = bisect_left(targets, power)
                sums[index] += log_prime
                reciprocal_sums[index] += log_prime / power
                counts[index] += 1
                power *= prime

        rows = []
        current_sum, current_reciprocal_sum, current_count = arb(0), arb(0), 0
        euler = arb.const_euler()
        for index, (x, thousandths) in enumerate(zip(targets, reported_millis, strict=True)):
            current_sum += sums[index]
            current_reciprocal_sum += reciprocal_sums[index]
            current_count += counts[index]
            result = x * current_reciprocal_sum - current_sum - x * arb(x).log() + (1 + euler) * x
            if not arb(2 * thousandths - 1) / 2000 < result < arb(2 * thousandths + 1) / 2000:
                raise ArithmeticError(f"G({x}) is outside the published rounding interval")
            rows.append({
                "x": x,
                "prime_powers_at_most_x": current_count,
                "published_thousandths": thousandths,
                "g_interval_dyadic": [dyadic(result.lower()), dyadic(result.upper())],
            })
    return {"precision_bits": bits, "rows": rows}


def run():
    """Reproduce the finite checks and the examples in III.5 and III.6."""
    example = {}
    for d0, d1, expected, positive_parts in ((1000, 10, Fraction(2405, 12), 190),
                                             (1000, 250, Fraction(625, 12), 0)):
        gamma = Fraction(5, 24)
        actual = gamma * (d0 + d1) - orbit_two_max_flow(5, d0, d1)
        # All odd powers from 5**3 onward give nonpositive terms in these examples.
        parts = max(0, Fraction(d0, 5) - d1)
        if actual != expected or parts != positive_parts or Fraction(d0, 5**3) > d1:
            raise ArithmeticError("the published orbit-two example does not replay")
        example[f"{d0},{d1}"] = {"F_exact": str(actual), "sum_of_positive_parts_exact": str(parts)}

    radius_example = {}
    for radius, expected in ((Fraction(1, 3), 0), (Fraction(3, 8), Fraction(-3, 8))):
        # In (III.125), the published example has A=9, prime=3, and J=1.
        slack = Fraction(27, 8) - Fraction(3, 8) - 9 * radius
        if slack != expected:
            raise ArithmeticError("the published radius example does not replay")
        radius_example[str(radius)] = str(slack)

    return {
        "schema": "NUMERIA.P93V40.RUBINSTEIN_REPLAY.1",
        "source": SOURCE,
        "scope": "Independent finite checks of the cited draft with rational and interval arithmetic.",
        "top_radius_III_92_and_93_rational_cases": check_top_radius(),
        "two_site_flow_III_97_rational_cases": check_orbit_two(),
        "five_adic_orbit_two_examples": example,
        "radius_example_III_125_slack_by_radius": radius_example,
        "smoothed_von_mangoldt_sum_I_2": check_prime_power_sums(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group()
    io.add_argument("--verify", type=Path)
    io.add_argument("--out", type=Path)
    parser.add_argument("--source-pdf", type=Path, help="also check the SHA-256 of a local v1 PDF")
    args = parser.parse_args()
    if args.source_pdf is not None and sha256(args.source_pdf.read_bytes()).hexdigest() != SOURCE["sha256"]:
        raise ValueError("the source PDF is not the pinned v1 draft")
    payload = (json.dumps(run(), indent=2) + "\n").encode("utf-8")
    if args.verify is not None:
        if args.verify.read_bytes() != payload:
            raise ValueError("the Rubinstein replay differs from the certificate")
        print(f"Verified {args.verify}")
    elif args.out is not None:
        args.out.write_bytes(payload)
    else:
        print(payload.decode("utf-8"), end="")


if __name__ == "__main__":
    main()
