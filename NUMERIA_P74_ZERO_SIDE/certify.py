"""Certify Suzuki's screw function on a closed rational interval in Arb."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from NUMERIA_ES8P259.suzuki import _dyadic
from .tails import HEIGHT, TailBounds
from .zeros import verified_ordinates


SCHEMA = "NUMERIA.P74.ZERO_SIDE.1"
MAX_DEPTH = 12


def rational(value):
    if type(value) not in (int, str, Fraction):
        raise ValueError("expected an exact rational as a string, integer or Fraction")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as err:
        raise ValueError("invalid exact rational") from err


def ball(point):
    return arb(point.numerator) / point.denominator


def secant_lower(left_value, right_value, width, terms):
    """Minimize the lower secant parabola on the entire closed cell."""
    left, right = left_value.lower(), right_value.lower()
    endpoint_min = min(left, right)
    curvature = terms * width**2
    slope = abs(right - left)
    if slope >= curvature:
        return endpoint_min
    return (endpoint_min - (curvature - slope)**2 / (4 * curvature)).lower()


class Spectrum:
    def __init__(self, ordinates):
        if not ordinates[-1] < HEIGHT:
            raise ArithmeticError("retained zeros must fall within the verified finite height")
        self.weights = [(gamma / 2, 4 / gamma**2) for gamma in ordinates]

    def value(self, t):
        time = ball(t)
        return sum((weight * (phase * time).sin()**2 for phase, weight in self.weights), arb(0))


def coverage(start, end, spectrum, tails, method, width, landmarks=()):
    count = len(spectrum.weights)
    evaluations = cells = attempts = deepest = 0
    total = None
    passed = start
    digest = sha256()

    def at(t):
        nonlocal evaluations
        evaluations += 1
        return spectrum.value(t)

    def cover(a, b, a_value, b_value, depth=0):
        nonlocal cells, attempts, deepest, passed, total
        attempts += 1
        cap = tails.bound(ball(b))
        margin = (secant_lower(a_value, b_value, ball(b - a), count) - cap).lower()
        cell = [str(a), str(b)]
        if margin > 0:
            cells += 1
            deepest = max(deepest, depth)
            passed = b
            if total is None or margin < total[0]:
                total = margin, cell
            record = [*cell, _dyadic(margin)]
            digest.update((json.dumps(record, separators=(",", ":")) + "\n").encode("ascii"))
            return None
        if method == "uniform" or depth >= MAX_DEPTH:
            return {
                "cell": cell,
                "cell_lower_dyadic": _dyadic(margin),
                "point_envelope_upper_dyadic": _dyadic((b_value - cap).upper()),
                "depth": depth,
            }

        mid = (a + b) / 2
        middle_value = at(mid)
        left_failure = cover(a, mid, a_value, middle_value, depth + 1)
        if left_failure is not None:
            return left_failure
        return cover(mid, b, middle_value, b_value, depth + 1)

    a, a_value = start, at(start)
    failure = None
    checkpoints = []
    targets = [value for value in landmarks if start < value < end] + [end]
    for target in targets:
        while a < target:
            b = min(target, a + width)
            b_value = at(b)
            failure = cover(a, b, a_value, b_value)
            if failure is not None:
                break
            a, a_value = b, b_value
        checkpoints.append({
            "target_exact": str(target),
            "status": "certified_positive" if failure is None else "inconclusive",
            "cells_through_target": cells if failure is None else None,
            "minimum_through_target_dyadic": _dyadic(total[0]) if failure is None else None,
        })
        if failure is not None:
            break

    endpoint_value = a_value if a == end else at(end)
    endpoint = (endpoint_value - tails.bound(ball(end))).upper()
    return {
        "status": "certified_positive" if failure is None else "inconclusive",
        "certified_through_exact": str(passed) if cells else None,
        "cell_count": cells,
        "attempted_cell_count": attempts,
        "endpoint_evaluation_count": evaluations,
        "maximum_subdivision_depth": deepest,
        "minimum": {"lower_dyadic": _dyadic(total[0]), "cell": total[1]} if total else None,
        "first_unresolved_cell": failure,
        "endpoint_envelope_upper_dyadic": _dyadic(endpoint),
        "checkpoints": checkpoints,
        "ordered_cell_bounds_sha256": digest.hexdigest(),
    }


def certify(start=18, end="60", zero_count=50, precision=192, tail_profile="graded",
            method="adaptive", mesh_width="1/64", checkpoints=()):
    start, end, width = rational(start), rational(end), rational(mesh_width)
    landmarks = tuple(rational(point) for point in checkpoints)
    if not 0 < start < end or not 0 < width <= 1:
        raise ValueError("expected 0 < start < end and 0 < mesh_width <= 1")
    if type(precision) is not int or precision < 96:
        raise ValueError("precision must be an integer of at least 96 bits")
    if type(zero_count) is not int or zero_count < 1:
        raise ValueError("zero_count must be a positive integer")
    if method not in ("adaptive", "uniform"):
        raise ValueError("method must be adaptive or uniform")
    if list(landmarks) != sorted(set(landmarks)) or not all(start < t < end for t in landmarks):
        raise ValueError("checkpoints must be increasing and strictly inside the interval")

    with ctx.workprec(precision):
        tails = TailBounds(tail_profile)
        ordinates, zero_metadata = verified_ordinates(zero_count)
        result = coverage(start, end, Spectrum(ordinates), tails, method, width, landmarks)
        return {
            "schema": SCHEMA,
            "start_exact": str(start),
            "end_exact": str(end),
            "zero_count": zero_count,
            "precision_bits": precision,
            "tail_profile": tail_profile,
            "method": method,
            "mesh_width_exact": str(width),
            "checkpoint_inputs_exact": [str(t) for t in landmarks],
            "overlaps_prime_side_to_10_pow_8": bool(ball(start).exp() < 10**8),
            "zeros": zero_metadata,
            "tail_bounds": tails.metadata(),
            **result,
        }


def verify_certificate(record):
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("unsupported P74 certificate schema")
    expected = certify(record.get("start_exact"), record.get("end_exact"),
                       record.get("zero_count"), record.get("precision_bits"),
                       record.get("tail_profile"), record.get("method"),
                       record.get("mesh_width_exact"), record.get("checkpoint_inputs_exact", ()))
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("certificate differs from rigorous replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    parser.add_argument("--start", default="18")
    parser.add_argument("--end", default="60")
    parser.add_argument("--zeros", type=int, default=50)
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--tail", choices=("coarse", "graded"), default="graded")
    parser.add_argument("--method", choices=("adaptive", "uniform"), default="adaptive")
    parser.add_argument("--mesh-width", default="1/64")
    args = parser.parse_args()
    if args.verify is not None:
        verify_certificate(json.loads(args.verify.read_text(encoding="utf-8")))
        print(f"Verified {args.verify}")
    else:
        result = certify(args.start, args.end, args.zeros, args.precision,
                         args.tail, args.method, args.mesh_width)
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
