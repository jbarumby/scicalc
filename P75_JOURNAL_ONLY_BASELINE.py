#!/usr/bin/env python3
"""Replay the P75 finite cover with only peer-reviewed tail estimates.

Run from the repository root on the audit branch (which contains P74):
  python P75_JOURNAL_ONLY_BASELINE.py --out /tmp/p75-journal.json

The KLN estimate counts beta > sigma0. For each target sigma, this script
uses the immediately preceding published table row sigma0 < sigma, so zeros
with beta >= sigma are included without assuming no zero lies on a grid line.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from flint import arb, ctx

from NUMERIA_ES8P259.suzuki import _dyadic
from NUMERIA_P74_ZERO_SIDE.certify import Spectrum, ball, rational
from NUMERIA_P74_ZERO_SIDE.tails import HEIGHT, count_tail_upper
from NUMERIA_P74_ZERO_SIDE.zeros import verified_ordinates

# Kadiri-Lumley-Ng, JMAA 465 (2018), Table 1, (4.71):
# N(sigma,T) <= A (log(kT))^(2 sigma) (log T)^(5-4 sigma)
#                    T^(8(1-sigma)/3) + B (log T)^2, with k=1.
KLN_ROWS = (
    (60, "2.177", "5.663"), (65, "2.963", "5.249"),
    (70, "3.983", "4.824"), (75, "5.277", "4.403"),
    (80, "6.918", "3.997"), (85, "8.975", "3.588"),
    (86, "9.441", "3.514"), (87, "9.926", "3.430"),
    (88, "10.431", "3.346"), (89, "10.955", "3.262"),
    (90, "11.499", "3.186"), (91, "12.063", "3.102"),
    (92, "12.646", "3.017"), (93, "13.250", "2.941"),
    (94, "13.872", "2.856"), (95, "14.513", "2.772"),
    (96, "15.173", "2.694"), (97, "15.850", "2.609"),
    (98, "16.544", "2.532"), (99, "17.253", "2.446"),
)
TARGETS = (80, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99)
SEGMENTS = 256
MAX_DEPTH = 12


def integrated_kln_density(sigma0, coefficient_a, coefficient_b):
    """Arb upper bound for 2*integral_H^inf N(sigma0,u) u^-3 du."""
    s = arb(sigma0) / 100
    log_h = arb(HEIGHT).log()
    alpha = arb(8) / 3 * (1 - s)
    shape = 6 - 2 * s
    decay = 2 - alpha
    if not decay > 0:
        raise ArithmeticError("the KLN logarithmic integral must converge")
    first = (decay * log_h).gamma_upper(shape) / decay**shape
    second = (log_h**2 + log_h + arb(1) / 2) / arb(HEIGHT)**2
    return (2 * arb(coefficient_a) * first + arb(coefficient_b) * second).upper()


class JournalOnlyTail:
    """Monotone cauda bound from peer-reviewed KLN, BW, and PT inputs."""
    def __init__(self):
        self.total_mass = count_tail_upper(HEIGHT).upper()
        self.theta = [arb(i) / (2 * SEGMENTS) for i in range(SEGMENTS + 1)]
        self.threshold_caps = []
        for target in TARGETS:
            # The source counts beta > sigma0; require sigma0 < target.
            sigma0, a, b = max(row for row in KLN_ROWS if row[0] < target)
            cap = min((self.total_mass / 2).upper(),
                      integrated_kln_density(sigma0, a, b))
            self.threshold_caps.append((arb(target) / 100, cap))

        self.mass_caps = []
        cap = (self.total_mass / 2).upper()
        index = 0
        for j in range(SEGMENTS):
            sigma = arb(1) / 2 + arb(j) / (2 * SEGMENTS)
            while (index < len(self.threshold_caps)
                   and self.threshold_caps[index][0] <= sigma):
                cap = min(cap, self.threshold_caps[index][1])
                index += 1
            self.mass_caps.append(cap)

    def bound(self, t):
        if not t >= 0:
            raise ValueError("t must be nonnegative")
        result = 2 * self.total_mass
        for j, cap in enumerate(self.mass_caps):
            left = self.theta[j] * t
            right = self.theta[j + 1] * t
            increment = 2 * ((left + right) / 2).sinh() * ((right - left) / 2).sinh()
            result += 4 * cap * increment.upper()
        coarse = (2 * self.total_mass * (t / 2).cosh()).upper()
        return min(result.upper(), coarse)


def local_curvature_upper(gammas, midpoint, width):
    """Upper-bound |S''| on a closed cell using S''' and Arb intervals."""
    second = (2 * sum(((g * midpoint).cos() for g in gammas), arb(0))).upper()
    third = (sum((g.upper() for g in gammas), arb(0)) * width.upper()).upper()
    return min(arb(2 * len(gammas)), (second + third).upper()).upper()


def lower_secant_local(left_value, right_value, width, curvature_upper):
    left, right = left_value.lower(), right_value.lower()
    endpoint_min = min(left, right)
    curvature = max(arb(0), curvature_upper.upper())
    if not curvature > 0:
        return endpoint_min
    c = (curvature * width.upper()**2 / 2).upper()
    d = abs(right - left).upper()
    if d >= c:
        return endpoint_min
    return (endpoint_min - ((c - d)**2 / (4 * c)).upper()).lower()


def certify(start, end, zero_count, precision, mesh):
    ordinates, zero_meta = verified_ordinates(zero_count)
    tails = JournalOnlyTail()
    spectrum = Spectrum(ordinates)
    count = len(ordinates)
    digest = __import__("hashlib").sha256()
    best = None
    accepted = attempts = 0
    cursor = start
    failure = None

    def value(t):
        return spectrum.value(t)

    def attempt(a, b, depth):
        nonlocal best, accepted, attempts, cursor
        attempts += 1
        midpoint = ball((a + b) / 2)
        width = ball(b - a)
        curvature = local_curvature_upper(ordinates, midpoint, width)
        secant = lower_secant_local(value(a), value(b), width, curvature)
        margin = (secant - tails.bound(ball(b))).lower()
        if margin > 0:
            accepted += 1
            cursor = b
            if best is None or margin < best[0]:
                best = margin, [str(a), str(b)]
            record = [str(a), str(b), _dyadic(margin), _dyadic(curvature)]
            digest.update((json.dumps(record, separators=(",", ":")) + "\n").encode("ascii"))
            return None
        if depth >= MAX_DEPTH:
            return {"cell": [str(a), str(b)], "lower_dyadic": _dyadic(margin), "depth": depth}
        mid = (a + b) / 2
        return attempt(a, mid, depth + 1) or attempt(mid, b, depth + 1)

    a = start
    while a < end:
        b = min(end, a + mesh)
        failure = attempt(a, b, 0)
        if failure is not None:
            break
        a = b

    return {
        "schema": "NUMERIA.P75.JOURNAL_ONLY_BASELINE.1",
        "tail_profile": "KLN-2018 strict-convention-safe + BW-2025 + PT-2021",
        "begin_t": str(start), "end_t": str(end),
        "zero_count": zero_count, "precision_bits": precision,
        "segments": SEGMENTS, "mesh_exact": str(mesh),
        "zero_provenance": zero_meta,
        "status": "conditional_finite_positive" if failure is None else "inconclusive",
        "certified_through_exact": str(cursor) if accepted else None,
        "cell_count": accepted, "attempted_cells": attempts,
        "minimum": {"lower_dyadic": _dyadic(best[0]), "cell": best[1]} if best else None,
        "first_unresolved_cell": failure,
        "ordered_cells_sha256": digest.hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", default="18")
    parser.add_argument("--end", default="79")
    parser.add_argument("--zeros", type=int, default=500)
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--mesh", default="1/64")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    start, end, mesh = map(rational, (args.start, args.end, args.mesh))
    if not (0 < start < end and 0 < mesh <= 1):
        parser.error("require 0 < start < end and 0 < mesh <= 1")
    with ctx.workprec(args.precision):
        result = certify(start, end, args.zeros, args.precision, mesh)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "status", "certified_through_exact", "cell_count",
        "attempted_cells", "minimum", "first_unresolved_cell")}, indent=2))


if __name__ == "__main__":
    main()
