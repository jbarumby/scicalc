"""Generate interval GL60 chunks for the reduced NUMERIA finite head."""

import argparse
import gzip
import json
from fractions import Fraction
from operator import index
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

from NUMERIA_ES8P249.legendre_tail import ARCHIVED_CUTOFF, HALF_WIDTH
from NUMERIA_ES8P249.prime_operator import P236_PRIME_POWERS


MODES = 259
CUTOFF = 309
SCHEMA = "NUMERIA.ES8P249.FINITE_HEAD_GL60_CHUNK.1"
SECTORS = ("even", "odd")


def _prime_terms():
    return tuple(
        (2 * arb(p).log() / arb(q).sqrt(), arb(q).log())
        for q, p in P236_PRIME_POWERS
    )


def _beta(prime_terms, pi):
    cutoff = arb(ARCHIVED_CUTOFF)
    mass = sum((amplitude for amplitude, _ in prime_terms), arb(0))
    return (cutoff / (2 * pi)).log() - 1 / cutoff - mass


def fixed_beta():
    """Enclose the original beta at the active Arb precision."""
    return _beta(_prime_terms(), arb.pi())


def spherical_bessel_j_sequence(x, max_degree):
    """Enclose all j_n(x) by downward recurrence from two native Arb seeds.

    The recurrence uses the identity j_(n-1)=(2n+1)j_n/x-j_(n+1).
    Both starting values are interval Bessel evaluations of the exact
    requested orders, so there is no truncation or Miller normalization.
    """
    max_degree = index(max_degree)
    if max_degree < 0:
        raise ValueError("max_degree must be nonnegative")
    x = arb(x)
    if not x.is_finite() or not x.lower() > 0:
        raise ValueError("the Bessel argument interval must be positive and finite")
    inverse_x = 1 / x
    scale = (arb.pi() * inverse_x / 2).sqrt()
    order = arb(2 * max_degree + 1) / 2
    values = [arb(0)] * (max_degree + 1)
    values[max_degree] = scale * x.bessel_j(order)
    following = scale * x.bessel_j(order + 1)
    for n in range(max_degree, 0, -1):
        values[n - 1] = (2 * n + 1) * inverse_x * values[n] - following
        following = values[n]
    if not values[0].overlaps(x.sin() * inverse_x):
        raise ArithmeticError("downward Bessel recurrence disagrees with j_0")
    return values


def _parameters(start, stop, cutoff, modes, bits):
    start, stop, cutoff, modes, bits = map(index, (start, stop, cutoff, modes, bits))
    if not 1 <= cutoff <= ARCHIVED_CUTOFF or not 0 <= start < stop <= 2 * cutoff:
        raise ValueError("require a nonempty panel range within the given cutoff")
    if not 1 <= modes <= 377 or bits < 128:
        raise ValueError("use 1 to 377 modes and at least 128 Arb matrix bits")
    return start, stop, cutoff, modes, bits


def _batch(start, stop, modes, bits, bessel_bits, roots, widths, norms, terms, pi, beta, vectors):
    even, odd, weighted_even, weighted_odd = [], [], [], []
    half_width, quarter = widths
    with ctx.workprec(bessel_bits):
        for panel in range(start, stop):
            center = arb(2 * panel + 1) / 4
            for root, weight in roots:
                t = center + quarter * root
                vals = spherical_bessel_j_sequence(half_width * t, 2 * modes - 1)
                prime = sum((a * (t * shift).cos() for a, shift in terms), arb(0))
                digamma = acb(arb(1) / 4, t / 2).digamma().real
                weight = quarter * weight * (digamma - pi.log() - prime - beta) / pi
                with ctx.workprec(bits):
                    corrected = [norm * value for norm, value in zip(norms, vals)]
                    evens, odds = corrected[0::2], corrected[1::2]
                    even.extend(evens)
                    odd.extend(odds)
                    weighted_even.extend(weight * value for value in evens)
                    weighted_odd.extend(weight * value for value in odds)

    rows = (stop - start) * 60
    outputs = []
    for parity, (u, wu) in enumerate(((even, weighted_even), (odd, weighted_odd))):
        U = arb_mat(rows, modes, u)
        WU = arb_mat(rows, modes, wu)
        weighted = WU if vectors is None else WU * vectors[parity]
        outputs.append(U.transpose() * weighted)
    return tuple(outputs)


def quadrature_chunk(
    start, stop, *, cutoff=CUTOFF, modes=MODES, bits=160,
    bessel_bits=640, threads=1, batch_panels=8, vectors=None, half_width=HALF_WIDTH,
):
    """Return the even and odd interval GL60 matrices for panels [start,stop).

    Each panel has width 1/2 and 60 nodes. The multiplier already subtracts
    the P237 beta at cutoff 504 and uses its four prime powers for every
    exact rational half width. The beta identity, the signed pole matrix,
    and the analytical quadrature remainder are added in the head certificate.
    The batch size bounds intermediate storage independently of the range.
    If given two vectors, return the enclosed GL60 actions as columns and
    form no dense Gram matrix. This evaluates U^T(WU x) on each batch.
    """
    start, stop, cutoff, modes, bits = _parameters(start, stop, cutoff, modes, bits)
    bessel_bits, threads, batch_panels = map(index, (bessel_bits, threads, batch_panels))
    half_width = Fraction(half_width)
    if half_width <= 0:
        raise ValueError("the half width must be positive")
    if bessel_bits < bits or threads < 1 or batch_panels < 1:
        raise ValueError("invalid Bessel precision, thread count, or batch size")
    old_threads = ctx.threads
    try:
        ctx.threads = threads
        with ctx.workprec(bessel_bits):
            roots = tuple(arb.legendre_p_root(60, k, weight=True) for k in range(60))
            if not all(roots[i][0].lower() > roots[i + 1][0].upper() for i in range(59)):
                raise ArithmeticError("Arb GL60 roots are not separated")
            if not all(weight.lower() > 0 for _, weight in roots):
                raise ArithmeticError("Arb GL60 weights are not positive")
            if not sum((weight for _, weight in roots), arb(0)).contains(2):
                raise ArithmeticError("Arb GL60 weight sum misses 2")
            width_ball = arb(half_width.numerator) / half_width.denominator
            widths = width_ball, arb(1) / 4
            norms = tuple(
                (-1) ** (n // 2) * (2 * width_ball * (2 * n + 1)).sqrt()
                for n in range(2 * modes)
            )
            terms = _prime_terms()
            pi = arb.pi()
            beta = _beta(terms, pi)

        with ctx.workprec(bits):
            columns = modes
            if vectors is not None:
                if len(vectors) != len(SECTORS) or any(len(vector) != modes for vector in vectors):
                    raise ValueError("expected one vector of the specified dimension per parity")
                vectors = tuple(arb_mat(modes, 1, list(vector)) for vector in vectors)
                if any(not entry.is_finite() for vector in vectors for entry in vector.entries()):
                    raise ValueError("expected finite coefficients in both vectors")
                columns = 1
            even, odd = arb_mat(modes, columns), arb_mat(modes, columns)
            for left in range(start, stop, batch_panels):
                right = min(left + batch_panels, stop)
                partial_even, partial_odd = _batch(
                    left, right, modes, bits, bessel_bits, roots, widths, norms, terms, pi, beta, vectors,
                )
                even += partial_even
                odd += partial_odd
            return even, odd
    finally:
        ctx.threads = old_threads


def _dyadic(endpoint):
    mantissa, exponent = endpoint.man_exp()
    return [str(mantissa), str(exponent)]


def chunk_report(start, stop, *, cutoff=CUTOFF, modes=MODES, bits=160, half_width=HALF_WIDTH, **options):
    """Serialize a computed chunk with exact, outward dyadic endpoints."""
    start, stop, cutoff, modes, bits = _parameters(start, stop, cutoff, modes, bits)
    half_width = Fraction(half_width)
    matrices = quadrature_chunk(
        start, stop, cutoff=cutoff, modes=modes, bits=bits, half_width=half_width, **options,
    )
    sectors = {}
    with ctx.workprec(bits):
        for sector, matrix in zip(SECTORS, matrices):
            entries = []
            for i in range(modes):
                for j in range(i, modes):
                    value = matrix[i, j].intersection(matrix[j, i])
                    if not value.is_finite():
                        raise ArithmeticError(f"nonfinite interval at {sector}[{i},{j}]")
                    entries.append([i, j, [_dyadic(value.lower()), _dyadic(value.upper())]])
            sectors[sector] = entries
    return {
        "schema": SCHEMA,
        "panel_range": [start, stop],
        "cutoff": cutoff,
        "beta_cutoff": ARCHIVED_CUTOFF,
        "half_width": str(half_width),
        "modes_per_parity": modes,
        "matrix_precision_bits": bits,
        "bessel_precision_bits": options.get("bessel_bits", 640),
        "GL60_nodes": 60 * (stop - start),
        "beta_identity_poles_and_analytic_remainder_added": False,
        "sectors": sectors,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--cutoff", type=int, default=CUTOFF)
    parser.add_argument("--modes", type=int, default=MODES)
    parser.add_argument("--half-width", type=Fraction, default=HALF_WIDTH)
    parser.add_argument("--bits", type=int, default=160)
    parser.add_argument("--bessel-bits", type=int, default=640)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--batch-panels", type=int, default=8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = chunk_report(
        args.start, args.stop, cutoff=args.cutoff, modes=args.modes, bits=args.bits, half_width=args.half_width,
        bessel_bits=args.bessel_bits, threads=args.threads, batch_panels=args.batch_panels,
    )
    args.out.write_bytes(gzip.compress(json.dumps(report, separators=(",", ":")).encode(), 3, mtime=0))


if __name__ == "__main__":
    main()
