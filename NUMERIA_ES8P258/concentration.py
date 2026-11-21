"""Directed GL60 and endpoint enclosures of Fourier–Legendre concentration."""

import argparse
from fractions import Fraction
import gzip
from hashlib import file_digest
import json
from math import isqrt
from operator import index
from pathlib import Path

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P251.gate_b.finite_head import _decode_sectors, _read_document, _require_header
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import _dyadic, spherical_bessel_j_sequence


HALF_WIDTH = Fraction("0.972955074527656652552676371721590")
CUTOFF = 504
MODES = 377
SECTORS = ("even", "odd")
CHUNK_SCHEMA = "NUMERIA.ES8P258.CONCENTRATION.GL60.CHUNK.1"
AGGREGATE_SCHEMA = "NUMERIA.ES8P258.CONCENTRATION.GL60.AGGREGATE.1"


def _geometry(half_width, cutoff, modes):
    half_width = Fraction(half_width)
    cutoff, modes = index(cutoff), index(modes)
    if not 0 < half_width <= 2 or cutoff < 1 or modes < 1:
        raise ValueError("invalid concentration geometry")
    return half_width, cutoff, modes


def _ceil_sqrt(value):
    floor = isqrt(value.numerator // value.denominator)
    return floor + (floor * floor < value)


def gl60_error_budget(sector, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES):
    """Bound the analytic error of the unweighted concentration quadrature."""
    half_width, cutoff, modes = _geometry(half_width, cutoff, modes)
    parity = SECTORS.index(sector)
    radius_squared = half_width**2 * ((cutoff + Fraction(1, 6))**2 + Fraction(1, 9))
    argument_upper = _ceil_sqrt(radius_squared)
    envelope = []
    factorial = 1
    for n in range(2 * modes - 1 + parity):
        factorial *= 2 * n + 1
        if n % 2 == parity:
            normalization = _ceil_sqrt(2 * half_width * (2 * n + 1))
            bessel = min(Fraction(2), Fraction(2 * argument_upper**n, factorial))
            envelope.append(normalization * bessel)

    # E3 quadrature bound 12/3^120, Jacobian 1/4, and 333/106 < pi.
    factor = Fraction(2 * cutoff * 106, 111 * 3**120)
    return {
        "factor": factor,
        "degree_envelopes": tuple(envelope),
        "row_sum_bound": factor * max(envelope) * sum(envelope),
    }


def concentration_chunk(
    start, stop, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES,
    bits=256, bessel_bits=1536, batch_panels=8,
):
    """Enclose both GL60 sums over half-unit panels [start, stop)."""
    half_width, cutoff, modes = _geometry(half_width, cutoff, modes)
    if not 0 <= start < stop <= 2 * cutoff or not 128 <= bits <= bessel_bits or batch_panels < 1:
        raise ValueError("invalid concentration panel range or precision")
    with ctx.workprec(bessel_bits):
        pi = arb.pi()
        width = arb(half_width.numerator) / half_width.denominator
        roots = tuple(arb.legendre_p_root(60, k, weight=True) for k in range(60))
        if not all(roots[k][0] > roots[k + 1][0] for k in range(59)):
            raise ArithmeticError("unseparated GL60 roots")
        if not all(weight > 0 for _, weight in roots):
            raise ArithmeticError("nonpositive GL60 weight")
        nodes = [(root / 4, (weight / (4 * pi)).sqrt()) for root, weight in roots]
        norms = [(-1)**(n // 2) * (2 * width * (2 * n + 1)).sqrt() for n in range(2 * modes)]

    with ctx.workprec(bits):
        totals = [arb_mat(modes, modes) for _ in SECTORS]
    for left in range(start, stop, batch_panels):
        right = min(left + batch_panels, stop)
        rows = [[], []]
        with ctx.workprec(bessel_bits):
            for panel in range(left, right):
                center = arb(2 * panel + 1) / 4
                for offset, weight in nodes:
                    values = spherical_bessel_j_sequence(width * (center + offset), 2 * modes - 1)
                    with ctx.workprec(bits):
                        row = [weight * norm * value for norm, value in zip(norms, values, strict=True)]
                    rows[0].extend(row[0::2])
                    rows[1].extend(row[1::2])
        with ctx.workprec(bits):
            for parity, row in enumerate(rows):
                block = arb_mat(60 * (right - left), modes, row)
                totals[parity] += block.transpose() * block
    return tuple(totals)


def _endpoint_concentration(panel, half_width, modes):
    result = tuple(arb_mat(modes, modes) for _ in SECTORS)
    if panel == 0:
        return result
    z = arb(half_width.numerator) * panel / (2 * half_width.denominator)
    pi = arb.pi()
    scale = (pi / (2 * z)).sqrt()
    values = [scale * z.bessel_j(arb(2 * n + 1) / 2) for n in range(2 * modes + 1)]
    derivatives = [n * values[n] / z - values[n + 1] for n in range(2 * modes)]
    norms = [(-1)**(n // 2) * arb(2 * n + 1).sqrt() for n in range(2 * modes)]

    diagonals = [(2 * z).si() - z.sin()**2 / z]
    for n in range(2 * modes - 1):
        boundary = z * (values[n]**2 + values[n + 1]**2)
        diagonals.append(((2 * n + 1) * diagonals[-1] - boundary) / (2 * n + 3))

    scale = 2 / pi
    for parity, matrix in enumerate(result):
        for i in range(modes):
            n = 2 * i + parity
            matrix[i, i] = scale * (2 * n + 1) * diagonals[n]
            for j in range(i + 1, modes):
                m = 2 * j + parity
                wronskian = values[m] * derivatives[n] - values[n] * derivatives[m]
                integral = z**2 * wronskian / (n * (n + 1) - m * (m + 1))
                value = scale * norms[n] * norms[m] * integral
                matrix[i, j] = matrix[j, i] = value
    if not all(value.is_finite() for matrix in result for value in matrix.entries()):
        raise ArithmeticError("nonfinite analytic concentration enclosure")
    return result


def analytic_concentration(
    start=0, stop=None, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES, bits=1536,
):
    """Enclose continuous G via the spherical Bessel ODE and sine integral."""
    half_width, cutoff, modes = _geometry(half_width, cutoff, modes)
    stop = 2 * cutoff if stop is None else index(stop)
    if not 0 <= start < stop <= 2 * cutoff or bits < 128:
        raise ValueError("invalid endpoint range or precision")
    with ctx.workprec(bits):
        lower = _endpoint_concentration(start, half_width, modes)
        upper = _endpoint_concentration(stop, half_width, modes)
        return tuple(a - b for a, b in zip(upper, lower, strict=True))


def _header(half_width, cutoff, modes):
    return {
        "cutoff": cutoff,
        "half_width": str(half_width),
        "modes_per_parity": modes,
        "analytic_remainder_added": False,
    }


def matrix_rows(matrix):
    for i in range(matrix.nrows()):
        for j in range(i, matrix.ncols()):
            value = matrix[i, j].intersection(matrix[j, i])
            if not value.is_finite():
                raise ArithmeticError("nonfinite concentration interval")
            yield [i, j, [_dyadic(value.lower()), _dyadic(value.upper())]]


def _write(path, header, matrices, bits):
    with ctx.workprec(bits):
        header["sectors"] = dict(zip(SECTORS, (list(matrix_rows(m)) for m in matrices), strict=True))
    data = json.dumps(header, separators=(",", ":"), allow_nan=False).encode()
    Path(path).write_bytes(gzip.compress(data, 3, mtime=0))


def write_chunk(path, start, stop, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES,
                bits=256, bessel_bits=1536):
    half_width, cutoff, modes = _geometry(half_width, cutoff, modes)
    blocks = concentration_chunk(
        start, stop, half_width=half_width, cutoff=cutoff, modes=modes, bits=bits, bessel_bits=bessel_bits,
    )
    header = _header(half_width, cutoff, modes) | {
        "schema": CHUNK_SCHEMA,
        "panel_range": [start, stop],
        "GL60_nodes": 60 * (stop - start),
        "matrix_precision_bits": bits,
        "bessel_precision_bits": bessel_bits,
    }
    _write(path, header, blocks, bits)


def aggregate_chunks(paths, path, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES):
    half_width, cutoff, modes = _geometry(half_width, cutoff, modes)
    boundary = 0
    totals = [arb_mat(modes, modes) for _ in SECTORS]
    records = []
    # Caller may list shards in any order; the content determines coverage.
    for shard in paths:
        doc, digest = _read_document(shard)
        _require_header(doc, _header(half_width, cutoff, modes) | {"schema": CHUNK_SCHEMA})
        start, stop = map(index, doc["panel_range"])
        if not 0 <= start < stop <= 2 * cutoff or doc["GL60_nodes"] != 60 * (stop - start):
            raise ValueError("invalid concentration chunk")
        records.append((start, stop, shard, digest))
    with ctx.workprec(384):
        manifest = []
        for start, stop, shard, digest in sorted(records):
            if start != boundary:
                raise ValueError("gap or overlap in concentration chunks")
            doc, actual_hash = _read_document(shard)
            if actual_hash != digest:
                raise ValueError("concentration chunk changed while loading")
            for i, block in enumerate(_decode_sectors(doc, modes)):
                totals[i] += block
            manifest.append({"panel_range": [start, stop], "sha256": digest})
            boundary = stop
    if boundary != 2 * cutoff:
        raise ValueError("incomplete concentration chunks")
    header = _header(half_width, cutoff, modes) | {
        "schema": AGGREGATE_SCHEMA,
        "panel_range": [0, 2 * cutoff],
        "GL60_nodes": 120 * cutoff,
        "chunks": manifest,
    }
    _write(path, header, totals, 384)


def load_authenticated(path, sha256, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES):
    half_width, cutoff, modes = _geometry(half_width, cutoff, modes)
    with Path(path).open("rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != sha256:
            raise ValueError("concentration aggregate SHA-256 mismatch")
    doc, digest = _read_document(path)
    if digest != sha256:
        raise ValueError("concentration aggregate changed while loading")
    _require_header(doc, _header(half_width, cutoff, modes) | {
        "schema": AGGREGATE_SCHEMA, "panel_range": [0, 2 * cutoff], "GL60_nodes": 120 * cutoff,
    })
    with ctx.workprec(384):
        return _decode_sectors(doc, modes)


def compare_concentration(quadrature, analytic, *, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES):
    """Require every endpoint interval to fit the GL60 enclosure and its error."""
    results = {}
    with ctx.workprec(384):
        for sector, sampled, continuous in zip(SECTORS, quadrature, analytic, strict=True):
            budget = gl60_error_budget(sector, half_width=half_width, cutoff=cutoff, modes=modes)
            factor = budget["factor"]
            scale = arb(factor.numerator) / factor.denominator
            envelope = [arb(q.numerator) / q.denominator for q in budget["degree_envelopes"]]
            max_radius = arb(0)
            for i in range(modes):
                for j in range(i, modes):
                    bound = (scale * envelope[i] * envelope[j]).upper()
                    exact = continuous[i, j]
                    outer = sampled[i, j] + arb(0, bound)
                    if not exact.is_finite() or not outer.is_finite() or not outer.contains(exact):
                        raise ArithmeticError(f"endpoint check failed at {sector}[{i},{j}]")
                    max_radius = max(max_radius, sampled[i, j].rad().upper())
            results[sector] = {
                "endpoint_enclosures_contained": modes * (modes + 1) // 2,
                "max_GL60_roundoff_radius_dyadic": _dyadic(max_radius),
            }
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--chunks", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.chunks:
        aggregate_chunks(args.chunks, args.out)
    elif args.start is not None and args.stop is not None:
        write_chunk(args.out, args.start, args.stop)
    else:
        parser.error("pass either --chunks or --start and --stop")


if __name__ == "__main__":
    main()
