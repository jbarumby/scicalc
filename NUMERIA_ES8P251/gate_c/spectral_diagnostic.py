"""Directed Rayleigh bounds for the authenticated Gate B continuous heads."""

import argparse
import json
from operator import index
from pathlib import Path

from flint import arb, arb_mat, ctx

from ..gate_b.finite_head import CUTOFF, MODES, SECTORS, load_aggregate
from ..gate_b.finite_head_certify import assemble_finite_head
from .bridge_audit import GATE_B_CERT_SHA, GATE_B_DATA, GATE_B_MATRIX_SHA, authenticated_input


DEFAULT_AGGREGATE = GATE_B_DATA / "P251_GATE_B_GL60.json.gz"
DEFAULT_CERTIFICATE = GATE_B_DATA / "P251_GATE_B_CERT.json"
_SEED_MODULUS = 13
_CANDIDATE_BITS = 80
_INVERSE_STEPS = 3
_POWER_STEPS = 32


def _point(value):
    mantissa, exponent = value.man_exp()
    return [str(mantissa), int(exponent)]


def _unwire_point(pair):
    mantissa, exponent = map(int, pair)
    return arb((mantissa, exponent))


def _rounded_column(vector):
    # The resulting points are exactly dyadic. Their quality only affects
    # the sharpness of the final variational bounds.
    with ctx.workprec(_CANDIDATE_BITS):
        return arb_mat(vector.nrows(), 1, [(value + 0).mid() for value in vector.entries()])


def _trial_vectors(matrix):
    center = matrix.mid()
    dimension = matrix.nrows()
    seed = arb_mat(
        dimension, 1,
        [arb((i * i + 2 * i + 5) % _SEED_MODULUS - 6) for i in range(dimension)],
    )
    low, high = seed, seed
    for _ in range(_INVERSE_STEPS):
        low = center.solve(low).mid()
    for _ in range(_POWER_STEPS):
        high = (center * high).mid()
    return _rounded_column(low), _rounded_column(high)


def _rayleigh(matrix, vector):
    if not all(value.is_exact() and value.is_finite() for value in vector.entries()):
        raise ValueError("the Rayleigh vector must have finite exact dyadic coordinates")
    squared_norm = (vector.transpose() * vector)[0, 0]
    if not squared_norm > 0:
        raise ValueError("the Rayleigh vector must be nonzero")
    quotient = (vector.transpose() * matrix * vector)[0, 0] / squared_norm
    if not quotient.is_finite():
        raise ArithmeticError("the Rayleigh quotient is nonfinite")
    return quotient


def _read_positive_certificate(path):
    citation = json.loads(authenticated_input(Path(path), GATE_B_CERT_SHA))
    if (
        citation["input_GL60_aggregate_sha256"] != GATE_B_MATRIX_SHA
        or citation["both_finite_heads_positive"] is not True
    ):
        raise ValueError("Gate B certificate does not bind both positive heads to this aggregate")
    lower_bounds = {}
    for sector in SECTORS:
        record = citation["sectors"][sector]
        if record["certified_positive"] is not True:
            raise ValueError(f"Gate B does not certify positivity in the {sector} sector")
        bound = _unwire_point(record["margin_lower_dyadic"])
        if not bound > 0:
            raise ValueError(f"Gate B has no positive margin in the {sector} sector")
        lower_bounds[sector] = bound
    return lower_bounds


def _sector_diagnostic(quadrature, sector, minimum_lower):
    matrix, metadata = assemble_finite_head(quadrature, sector, cutoff=CUTOFF)
    remainder = metadata["gl60_error_budget"]["row_sum_bound"]
    error = (arb(remainder.numerator) / remainder.denominator).upper()
    low, high = _trial_vectors(matrix)
    low_rayleigh = _rayleigh(matrix, low) + arb(0, error)
    high_rayleigh = _rayleigh(matrix, high) + arb(0, error)

    # Rayleigh's variational principle applies to the exact symmetric
    # continuous head. The interval matrix covers quadrature/assembly error.
    # The separately budgeted GL60 remainder shifts each quotient by at most
    # its operator-norm bound, and shifts the full norm by the same amount.
    minimum_upper = low_rayleigh.upper()
    maximum_lower = high_rayleigh.lower()
    row_bounds = [
        sum((entry.abs_upper() for entry in row), arb(0))
        for row in matrix.tolist()
    ]
    maximum_upper = (max(row.upper() for row in row_bounds) + error).upper()
    if not 0 < minimum_lower <= minimum_upper <= maximum_upper:
        raise ArithmeticError("inconsistent minimum eigenvalue bounds")
    if not 0 < maximum_lower <= maximum_upper:
        raise ArithmeticError("inconsistent maximum eigenvalue bounds")
    condition_lower = (maximum_lower / minimum_upper).lower()

    return {
        "lambda_min_lower_from_gate_b_dyadic": _point(minimum_lower),
        "lambda_min_upper_dyadic": _point(minimum_upper),
        "lambda_max_lower_dyadic": _point(maximum_lower),
        "lambda_max_upper_dyadic": _point(maximum_upper),
        "condition_number_lower_dyadic": _point(condition_lower),
        "GL60_analytic_operator_norm_upper_dyadic": _point(error),
        "low_trial_rayleigh_interval_dyadic": [_point(low_rayleigh.lower()), _point(low_rayleigh.upper())],
        "high_trial_rayleigh_interval_dyadic": [_point(high_rayleigh.lower()), _point(high_rayleigh.upper())],
        "low_trial_vector_dyadic": [_point(value) for value in low.entries()],
        "high_trial_vector_dyadic": [_point(value) for value in high.entries()],
    }


def run_diagnostic(
    matrix_path=DEFAULT_AGGREGATE, certificate_path=DEFAULT_CERTIFICATE, *, precision_bits=384,
):
    """Return reproducible, outward spectral-scale bounds for both Gate B heads."""
    precision_bits = index(precision_bits)
    if precision_bits < 128:
        raise ValueError("use at least 128 Arb bits")
    with ctx.workprec(precision_bits):
        bounds = _read_positive_certificate(certificate_path)
        matrices, _, digest = load_aggregate(Path(matrix_path))
        if digest != GATE_B_MATRIX_SHA:
            raise ValueError("Gate B quadrature aggregate SHA-256 mismatch")
        sectors = {
            sector: _sector_diagnostic(matrix, sector, bounds[sector])
            for sector, matrix in zip(SECTORS, matrices, strict=True)
        }
    return {
        "schema": "NUMERIA.ES8P251.GATE_C.SPECTRAL_DIAGNOSTIC.1",
        "input_gate_b_aggregate_sha256": GATE_B_MATRIX_SHA,
        "input_gate_b_positivity_certificate_sha256": GATE_B_CERT_SHA,
        "cutoff": CUTOFF,
        "modes_per_parity": MODES,
        "arb_precision_bits": precision_bits,
        "candidate_precision_bits": _CANDIDATE_BITS,
        "inverse_iteration_steps": _INVERSE_STEPS,
        "power_iteration_steps": _POWER_STEPS,
        "sectors": sectors,
        "scope": "continuous finite Legendre heads with the classical P237 GL60 budget",
        "Riemann_hypothesis_proved": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_AGGREGATE)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run_diagnostic(args.matrix, args.certificate, precision_bits=args.precision_bits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
