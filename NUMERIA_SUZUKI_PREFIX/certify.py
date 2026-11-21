"""Certify the prime-power prefixes of Suzuki's arithmetic inequality in Arb."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from NUMERIA_ES8P259.suzuki import DEFAULT_LIMIT as P259_LIMIT, _dyadic
from NUMERIA_SUZUKI_PREFIX.prime_powers import iter_prime_powers


SCHEMA = "NUMERIA.SUZUKI.PRIME_PREFIX.1"
DEFAULT_LIMIT = 10**8
DEFAULT_BLOCK = 2**20
DEFAULT_PRECISION = 192


def tangent_gap(a, b):
    """Enclose the minimum of 4*exp(t/2) - a*t + b over real t."""
    return b - 2 * a * ((a / 2).log() - 1)


def _quadratic_gap(a, b, root, logq):
    deficit = 4 * root - (a * logq - b)
    lower_bound = deficit - (a - 2 * root)**2 / root
    return deficit, lower_bound


def quadratic_gap(a, b, q):
    """Enclose the local deficit and its quadratic correction at q."""
    q = arb(q)
    return _quadratic_gap(a, b, q.sqrt(), q.log())


def _strict_minimum(previous, lower, q):
    if previous is None or lower < previous[0]:
        return lower, q
    return previous


def _minimum_record(minimum):
    if minimum is None:
        return None
    lower, q = minimum
    return {"prime_power": q, "lower_dyadic": _dyadic(lower)}


class _ScrewExtension:
    """Bound the archimedean correction above the last P259 breakpoint."""

    def __init__(self, start):
        quarter = arb(1) / 4
        self.shift = (arb.pi().log() - quarter.digamma()) / 2
        self.constant = arb(2).zeta(quarter) / 4 - 8
        if not self.shift > 0:
            raise ArithmeticError("could not enclose the positive slope shift")
        x = arb(start)
        self.tail_upper = (4 / (25 * (x*x - 1) * x.sqrt())).upper()
        self.minimum = None
        self.count = 0

    def add_prefix(self, a, b, q):
        correction = tangent_gap(a + self.shift, b + self.constant) - self.tail_upper
        self.minimum = _strict_minimum(self.minimum, correction.lower(), q)
        self.count += 1

    def record(self, limit):
        return {
            "from_log_of": P259_LIMIT,
            "through_log_of": limit,
            "active_prefix_count": self.count,
            "uniform_lerch_tail_upper_dyadic": _dyadic(self.tail_upper),
            "lower_bound_minimum": _minimum_record(self.minimum),
            "strictly_positive": bool(self.minimum[0] > 0),
        }


def certify(limit=DEFAULT_LIMIT, block_size=DEFAULT_BLOCK, precision=DEFAULT_PRECISION):
    """Check every prefix, retaining strict bounds and a P259-compatible digest."""
    if type(limit) is not int or limit < 2:
        raise ValueError("the certificate limit must be an integer of at least 2")
    if type(precision) is not int or precision < 64:
        raise ValueError("precision must be an integer of at least 64 bits")
    with ctx.workprec(precision):
        digest = sha256()
        overlap_digest = None
        count = new_count = 0
        minimum = strong_minimum = new_minimum = new_strong_minimum = None
        a, b = arb(0), arb(0)
        last_q = None
        screw = _ScrewExtension(P259_LIMIT) if limit > P259_LIMIT else None
        for q, p in iter_prime_powers(limit, block_size):
            if q > P259_LIMIT and overlap_digest is None:
                overlap_digest = digest.hexdigest()
                screw.add_prefix(a, b, last_q)

            point = arb(q)
            root, logq = point.sqrt(), point.log()
            weight = (logq if p == q else arb(p).log()) / root
            a += weight
            b += weight * logq

            gap = tangent_gap(a, b)
            lower = gap.lower()
            if not lower > 0:
                status = "strictly negative" if gap.upper() < 0 else "undecided"
                raise ArithmeticError(f"tangent gap at {q} is {status}: {gap}")
            minimum = _strict_minimum(minimum, lower, q)
            row = json.dumps([q, _dyadic(lower)], separators=(",", ":"))
            digest.update((row + "\n").encode("ascii"))

            _, strong = _quadratic_gap(a, b, root, logq)
            strong_lower = strong.lower()
            strong_minimum = _strict_minimum(strong_minimum, strong_lower, q)
            count += 1
            if q > P259_LIMIT:
                new_count += 1
                new_minimum = _strict_minimum(new_minimum, lower, q)
                new_strong_minimum = _strict_minimum(new_strong_minimum, strong_lower, q)
                screw.add_prefix(a, b, q)
            last_q = q

        if limit >= P259_LIMIT and overlap_digest is None:
            overlap_digest = digest.hexdigest()
        if screw is not None and screw.count == 0:
            screw.add_prefix(a, b, last_q)

        return {
            "schema": SCHEMA,
            "limit": limit,
            "precision_bits": precision,
            "prime_power_count": count,
            "tangent_minimum": _minimum_record(minimum),
            "quadratic_minimum": _minimum_record(strong_minimum),
            "quadratic_bound_positive_for_all_prefixes": bool(strong_minimum[0] > 0),
            "prefix_tangent_certificate_sha256": digest.hexdigest(),
            "p259_overlap_limit": P259_LIMIT if overlap_digest is not None else None,
            "p259_overlap_sha256": overlap_digest,
            "new_prime_power_count": new_count,
            "new_tangent_minimum": _minimum_record(new_minimum),
            "new_quadratic_minimum": _minimum_record(new_strong_minimum),
            "screw_extension": screw.record(limit) if screw is not None else None,
        }


def verify_certificate(record, block_size=DEFAULT_BLOCK):
    """Regenerate the evidence and reject any changed field or JSON type."""
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("unsupported certificate schema")
    expected = certify(record.get("limit"), block_size, record.get("precision_bits"))
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("certificate differs from rigorous replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path, help="write a new JSON certificate")
    io.add_argument("--verify", type=Path, help="repeat and compare a JSON certificate")
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--block-size", type=int, default=DEFAULT_BLOCK)
    parser.add_argument("--precision", type=int, default=DEFAULT_PRECISION)
    args = parser.parse_args()
    if args.verify is not None:
        verify_certificate(json.loads(args.verify.read_text(encoding="utf-8")), args.block_size)
        print(f"Verified {args.verify}")
    else:
        result = certify(args.limit, args.block_size, args.precision)
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
