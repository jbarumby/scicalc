"""Rebuild and certify the ES8-P251 Gate B finite head from interval GL60 chunks."""

import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re

from flint import arb, arb_mat, ctx

from .finite_head_certify import certify_finite_head
from .finite_head_quadrature import CUTOFF, MODES, SCHEMA, SECTORS, quadrature_chunk
from NUMERIA_ES8P249.interval_audit import certify_p236_tail_start
from NUMERIA_ES8P249.legendre_tail import ARCHIVED_CUTOFF, HALF_WIDTH, modal_tail_bounds


ASSEMBLY_BITS = 384
AGGREGATE_SCHEMA = "NUMERIA.ES8P251.GATE_B.DIRECT_GL60_HEAD.1"
CERTIFICATE_SCHEMA = "NUMERIA.ES8P251.GATE_B.FINITE_HEAD_CERTIFICATE.1"


def _read_document(path):
    raw = Path(path).read_bytes()
    return json.loads(gzip.decompress(raw)), sha256(raw).hexdigest()


def _require_header(document, expected):
    for field, value in expected.items():
        if field not in document or type(document[field]) is not type(value) or document[field] != value:
            raise ValueError(f"incompatible {field}")


def _parameters(cutoff, modes):
    return {
        "cutoff": cutoff,
        "beta_cutoff": ARCHIVED_CUTOFF,
        "half_width": str(HALF_WIDTH),
        "modes_per_parity": modes,
        "beta_identity_poles_and_analytic_remainder_added": False,
    }


def _decode_point(wire):
    if not isinstance(wire, list) or len(wire) != 2 or any(type(part) is not str for part in wire):
        raise ValueError("expected a dyadic encoded as two integers")
    mantissa, exponent = map(int, wire)
    return arb((mantissa, exponent))


def _decode_sectors(document, modes):
    sectors = document.get("sectors")
    if not isinstance(sectors, dict) or set(sectors) != set(SECTORS):
        raise ValueError("expected precisely the even and odd sectors")
    matrices = []
    for sector in SECTORS:
        entries = sectors[sector]
        if not isinstance(entries, list) or len(entries) != modes * (modes + 1) // 2:
            raise ValueError(f"incomplete upper triangle for {sector}")
        matrix = arb_mat(modes, modes)
        pairs = ((i, j) for i in range(modes) for j in range(i, modes))
        for (i, j), entry in zip(pairs, entries, strict=True):
            if (not isinstance(entry, list) or len(entry) != 3
                    or any(type(coordinate) is not int for coordinate in entry[:2]) or entry[:2] != [i, j]):
                raise ValueError(f"unexpected triangular entry for {sector} at {(i, j)}")
            bounds = entry[2]
            if not isinstance(bounds, list) or len(bounds) != 2:
                raise ValueError("expected both outward interval endpoints")
            lower, upper = map(_decode_point, bounds)
            if not lower <= upper:
                raise ValueError("inverted interval")
            matrix[i, j] = matrix[j, i] = lower.union(upper)
        matrices.append(matrix)
    return tuple(matrices)


def _check_coverage(records, cutoff):
    if not isinstance(records, list):
        raise ValueError("expected an ordered GL60 chunk manifest")
    end = 0
    for record in records:
        bounds = record.get("panel_range") if isinstance(record, dict) else None
        if (not isinstance(bounds, list) or len(bounds) != 2
                or any(type(endpoint) is not int for endpoint in bounds)):
            raise ValueError("invalid panel range in GL60 manifest")
        digest = record.get("sha256")
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("invalid GL60 chunk SHA-256")
        start, stop = bounds
        if start != end or not start < stop <= 2 * cutoff:
            raise ValueError("panel coverage has a gap, duplicate, overlap or invalid endpoint")
        end = stop
    if end != 2 * cutoff:
        raise ValueError("incomplete panel coverage")


def aggregate_chunks(paths, *, cutoff=CUTOFF, modes=MODES):
    """Sum complete, nonoverlapping GL60 chunks and bind the compressed bytes."""
    chunks = []
    for path in paths:
        document, digest = _read_document(path)
        _require_header(document, _parameters(cutoff, modes) | {"schema": SCHEMA})
        bounds = document.get("panel_range")
        if (not isinstance(bounds, list) or len(bounds) != 2
                or any(type(endpoint) is not int for endpoint in bounds)):
            raise ValueError("invalid panel range")
        start, stop = bounds
        _require_header(document, {"GL60_nodes": 60 * (stop - start)})
        matrix_bits = document.get("matrix_precision_bits")
        bessel_bits = document.get("bessel_precision_bits")
        if type(matrix_bits) is not int or type(bessel_bits) is not int or not bessel_bits >= matrix_bits >= 128:
            raise ValueError("invalid precision metadata")
        record = {
            "panel_range": bounds,
            "sha256": digest,
            "matrix_precision_bits": matrix_bits,
            "bessel_precision_bits": bessel_bits,
        }
        chunks.append((record, Path(path)))

    chunks.sort(key=lambda chunk: chunk[0]["panel_range"])
    records = [record for record, _ in chunks]
    _check_coverage(records, cutoff)
    with ctx.workprec(ASSEMBLY_BITS):
        even, odd = arb_mat(modes, modes), arb_mat(modes, modes)
        for record, path in chunks:
            document, digest = _read_document(path)
            if digest != record["sha256"]:
                raise ValueError("chunk changed during aggregation")
            partial_even, partial_odd = _decode_sectors(document, modes)
            even += partial_even
            odd += partial_odd
    return (even, odd), records


def _point(value):
    mantissa, exponent = value.man_exp()
    return [str(mantissa), str(exponent)]


def write_aggregate(path, matrices, records, *, cutoff=CUTOFF, modes=MODES):
    """Save the unaugmented sum for independent positivity replay."""
    sectors = {}
    with ctx.workprec(ASSEMBLY_BITS):
        for sector, matrix in zip(SECTORS, matrices, strict=True):
            sectors[sector] = [
                [i, j, [_point(matrix[i, j].lower()), _point(matrix[i, j].upper())]]
                for i in range(modes) for j in range(i, modes)
            ]
    document = _parameters(cutoff, modes) | {
        "schema": AGGREGATE_SCHEMA,
        "panel_range": [0, 2 * cutoff],
        "GL60_nodes": 120 * cutoff,
        "summation_precision_bits": ASSEMBLY_BITS,
        "chunks": records,
        "sectors": sectors,
    }
    raw = json.dumps(document, separators=(",", ":")).encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(raw, compresslevel=3, mtime=0))


def _dyadic_fraction(pair):
    mantissa, exponent = map(int, pair)
    return Fraction(mantissa) * Fraction(2)**exponent


def _schur_transfer(sectors):
    """Transfer a positive finite margin to the model using exact budgets."""
    if not all(report["certified_positive"] for report in sectors.values()):
        return {"certified_for_reduced_model": False}
    minimum_head = min(_dyadic_fraction(report["margin_lower_dyadic"]) for report in sectors.values())
    beta_lower = min(
        _dyadic_fraction(report["assembly"]["beta_endpoints_dyadic"][0])
        for report in sectors.values()
    )
    budgets = [modal_tail_bounds(CUTOFF, MODES, parity) for parity in (0, 1)]
    coupling_squared = max(pair[0] for pair in budgets)
    deviation = max(pair[1] for pair in budgets)
    # Apply Young with coefficient m/2 on the head and 2||B||²/m on the tail.
    floor = min(minimum_head / 2, beta_lower - deviation - 2 * coupling_squared / minimum_head)
    return {
        "certified_for_reduced_model": floor > 0,
        "model_floor_lower_exact": str(floor),
        "relies_on_P239_classical_modal_envelopes": True,
    }


def load_aggregate(path):
    """Load the complete Gate B aggregate and its compressed digest."""
    document, digest = _read_document(path)
    _require_header(document, _parameters(CUTOFF, MODES) | {
        "schema": AGGREGATE_SCHEMA,
        "panel_range": [0, 2 * CUTOFF],
        "GL60_nodes": 120 * CUTOFF,
        "summation_precision_bits": ASSEMBLY_BITS,
    })
    _check_coverage(document["chunks"], CUTOFF)
    with ctx.workprec(ASSEMBLY_BITS):
        matrices = _decode_sectors(document, MODES)
    return matrices, document["chunks"], digest


def certify_aggregate(path, *, precision_bits=1024):
    """Audit all 518 center pivots, Weyl perturbations and the model transfer."""
    matrices, chunks, digest = load_aggregate(path)
    sectors = {
        sector: certify_finite_head(matrix, sector, cutoff=CUTOFF, precision_bits=precision_bits)
        for sector, matrix in zip(SECTORS, matrices, strict=True)
    }
    tail = certify_p236_tail_start()
    if tail["certified_range"][0] > CUTOFF:
        raise ArithmeticError("the signed high-frequency certificate does not reach this cutoff")
    return {
        "schema": CERTIFICATE_SCHEMA,
        "input_GL60_aggregate_sha256": digest,
        "input_GL60_chunks": chunks,
        "L": str(HALF_WIDTH),
        "cutoff": CUTOFF,
        "beta_reference_cutoff": ARCHIVED_CUTOFF,
        "panels_per_sector": 2 * CUTOFF,
        "GL60_nodes_per_sector": 120 * CUTOFF,
        "modes_per_sector": MODES,
        "sectors": sectors,
        "both_finite_heads_positive": all(report["certified_positive"] for report in sectors.values()),
        "signed_high_frequency_tail": tail,
        "schur_transfer": _schur_transfer(sectors),
        "universal_Weil_positivity_certified": False,
        "Riemann_hypothesis_proved": False,
    }


def audit_matrix_free(path):
    """Compare full matrix-free GL60 actions to the independently summed matrix."""
    matrices, _, digest = load_aggregate(path)
    numerators = [[((i * i + 5 * i + 3 * parity) % 19) - 9 for i in range(MODES)]
                  for parity in (0, 1)]
    vectors = [[arb((numerator, -3)) for numerator in coefficients] for coefficients in numerators]
    actions = quadrature_chunk(0, 2 * CUTOFF, bits=256, vectors=vectors)
    sectors = {}
    with ctx.workprec(ASSEMBLY_BITS):
        for sector, matrix, vector, action in zip(SECTORS, matrices, vectors, actions, strict=True):
            reference = matrix * arb_mat(MODES, 1, vector)
            errors = []
            for direct, assembled in zip(action.entries(), reference.entries(), strict=True):
                difference = direct - assembled
                if not difference.contains(0):
                    raise ArithmeticError(f"matrix-free action disagrees with the {sector} GL60 aggregate")
                errors.append(difference.abs_upper())
            sectors[sector] = {"overlapping_coordinates": len(errors),
                               "difference_abs_upper_dyadic": _point(max(errors))}
    return {
        "schema": "NUMERIA.ES8P251.GATE_B.MATRIX_FREE_GL60_ACTION.1",
        "input_GL60_aggregate_sha256": digest,
        "panel_range": [0, 2 * CUTOFF],
        "vector_numerators_over_8": dict(zip(SECTORS, numerators, strict=True)),
        "sectors": sectors,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunks", type=Path, nargs="+", help="combine freshly generated GL60 chunks")
    parser.add_argument("--matrix", type=Path, required=True, help="aggregate GL60 matrix file")
    parser.add_argument("--out", type=Path, required=True, help="certificate JSON")
    parser.add_argument("--precision-bits", type=int, default=1024)
    parser.add_argument("--action-audit", type=Path, help="also replay a matrix-free action over all panels")
    args = parser.parse_args()
    if args.chunks is not None:
        matrices, records = aggregate_chunks(args.chunks)
        write_aggregate(args.matrix, matrices, records)
    report = certify_aggregate(args.matrix, precision_bits=args.precision_bits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for sector in SECTORS:
        print(sector, report["sectors"][sector]["status"], flush=True)
    if args.action_audit is not None:
        action = audit_matrix_free(args.matrix)
        args.action_audit.parent.mkdir(parents=True, exist_ok=True)
        args.action_audit.write_text(json.dumps(action, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
