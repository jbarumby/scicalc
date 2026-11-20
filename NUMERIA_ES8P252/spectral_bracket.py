"""Bracket the two P252 minimum eigenvalues by rigorous shifted sign tests."""

import argparse
import json
from pathlib import Path

from flint import arb, ctx

from NUMERIA_ES8P251.gate_b.finite_head import load_aggregate
from NUMERIA_ES8P251.gate_b.finite_head_certify import (
    _max_radius_row_sum,
    _midpoint_ldl_pivots,
    assemble_finite_head,
)
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import MODES, SECTORS
from NUMERIA_ES8P251.gate_c.spectral_diagnostic import _point, _rayleigh, _trial_vectors

from .window_audit import HALF_WIDTH


AGGREGATE_SHA256 = "85cccdafd673222792b59ed1e1c0d25f3c39586fbac808958c4af9619c65bbe7"
AGGREGATE = Path(__file__).resolve().parent / "data/P252_GL60.json.gz"
PRECISION_BITS = 1024


def check_shifted_center(center, lower, perturbation):
    """Certify a proposed lower bound with an LDL of the shifted center."""
    shifted = center.mid()
    shift = lower + perturbation
    for i in range(shifted.nrows()):
        shifted[i, i] -= shift
    pivots, failed = _midpoint_ldl_pivots(shifted)
    if failed is not None:
        raise ArithmeticError("shifted interval LDL did not certify every pivot positive")
    return pivots


def bracket_smallest(matrix, trial, gl60_error):
    """Return exact dyadic endpoints using a Rayleigh witness and shifted LDL.

    The supplied symmetric matrix encloses the quadrature assembly. The
    continuous GL60 error is supplied as a separate operator norm bound.
    """
    if gl60_error < 0:
        raise ValueError("the GL60 operator norm bound must be nonnegative")
    analytic = (arb(gl60_error.numerator) / gl60_error.denominator).upper()
    perturbation = (_max_radius_row_sum(matrix) + analytic).upper()
    upper = (_rayleigh(matrix, trial) + analytic).upper()

    # The trial need not be an eigenvector. Try wider brackets if its
    # Rayleigh value leaves too little room for the proposed lower shift.
    step = max(perturbation, arb((1, -(ctx.prec // 2))))
    for factor in (4, 16, 64, 256):
        lower = (upper - factor * step).lower()
        try:
            pivots = check_shifted_center(matrix, lower, perturbation)
        except ArithmeticError:
            continue
        return {
            "lambda_min_strict_lower_dyadic": _point(lower),
            "lambda_min_upper_dyadic": _point(upper),
            "bracket_width_upper_dyadic": _point((upper - lower).upper()),
            "quadrature_and_analytic_norm_error_upper_dyadic": _point(perturbation),
            "shifted_LDL_positive_pivots": len(pivots),
            "shifted_LDL_minimum_pivot_lower_dyadic": _point(min(pivot.lower() for pivot in pivots)),
        }
    raise ArithmeticError("the Rayleigh trial did not yield a certifiable bracket at the attempted shifts")


def run_audit(path=AGGREGATE):
    matrices, _, digest = load_aggregate(path, modes=MODES, half_width=HALF_WIDTH)
    if digest != AGGREGATE_SHA256:
        raise ArithmeticError("the P252 aggregate failed SHA-256 authentication")

    sectors = {}
    with ctx.workprec(PRECISION_BITS):
        for sector, matrix in zip(SECTORS, matrices, strict=True):
            head, metadata = assemble_finite_head(matrix, sector, half_width=HALF_WIDTH)
            trial, _ = _trial_vectors(head)
            error = metadata["gl60_error_budget"]["row_sum_bound"]
            result = bracket_smallest(head, trial, error)
            result["trial_vector_exact_dyadic"] = [_point(value) for value in trial.entries()]
            sectors[sector] = result
    return {
        "schema": "NUMERIA.ES8P252.SHIFTED_SPECTRAL_BRACKET.1",
        "input_GL60_aggregate_sha256": digest,
        "half_width_exact": str(HALF_WIDTH),
        "modes_per_sector": MODES,
        "precision_bits": PRECISION_BITS,
        "sectors": sectors,
        "scope": "continuous finite Legendre heads under the archived analytic GL60 bound",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=AGGREGATE)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run_audit(args.matrix)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
