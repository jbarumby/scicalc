"""Certify negative Rayleigh witnesses for the early four-prime auxiliary trial."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P249.legendre_tail import MAX_FOUR_PRIME_MODEL_HALF_WIDTH
from NUMERIA_ES8P251.gate_b.finite_head import load_aggregate
from NUMERIA_ES8P251.gate_b.finite_head_certify import _midpoint_ldl_pivots, assemble_finite_head
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import CUTOFF, MODES, SECTORS
from NUMERIA_ES8P251.gate_c.spectral_diagnostic import _point, _rayleigh


HALF_WIDTH = MAX_FOUR_PRIME_MODEL_HALF_WIDTH
EXPECTED_SHA256 = "f24ebe6d3bfd0e5c0d751c31267251c6a7b28dd0bab1dd14917487b12d182f4c"


def certify_negative_heads(path, *, expected_digest, cutoff=CUTOFF, half_width=HALF_WIDTH, modes=MODES):
    quadratures, _, digest = load_aggregate(path, half_width=half_width, cutoff=cutoff, modes=modes)
    if digest != expected_digest:
        raise ArithmeticError("low-cutoff trial matrix SHA-256 mismatch")
    results = {}
    with ctx.workprec(1024):
        for sector, quadrature in zip(SECTORS, quadratures, strict=True):
            head, assembly = assemble_finite_head(quadrature, sector, half_width=half_width, cutoff=cutoff)
            pivots, failed = _midpoint_ldl_pivots(head.mid())
            if failed is None or not failed < 0:
                raise ArithmeticError(f"could not locate a negative pivot for {sector}")
            k = len(pivots)
            leading = arb_mat(k, k, [head[i, j].mid() for i in range(k) for j in range(k)])
            column = arb_mat(k, 1, [-head[i, k].mid() for i in range(k)])
            trial = [entry.mid() for entry in leading.solve(column).entries()] if k else []
            trial.append(arb(1))
            vector = arb_mat(head.nrows(), 1, trial + [arb(0)] * (head.nrows() - len(trial)))
            error = assembly["gl60_error_budget"]["row_sum_bound"]
            upper = (_rayleigh(head, vector) + arb(error.numerator) / error.denominator).upper()
            if not upper < 0:
                raise ArithmeticError(f"negative Rayleigh quotient not certified for {sector}")
            results[sector] = {
                "vector_nonpadding_coordinates_exact_dyadic": [_point(value) for value in trial],
                "vector_dimension": head.nrows(),
                "continuous_Rayleigh_upper_dyadic": _point(upper),
            }
    return {
        "schema": "NUMERIA.ES8P253.LOW_CUTOFF_WITNESSES.1",
        "scope": "four-prime truncated auxiliary model",
        "half_width": str(half_width),
        "cutoff": cutoff,
        "modes_per_sector": modes,
        "input_GL60_aggregate_sha256": digest,
        "sectors": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sha256", default=EXPECTED_SHA256)
    parser.add_argument("--cutoff", type=int, default=CUTOFF)
    parser.add_argument("--half-width", type=Fraction, default=HALF_WIDTH)
    parser.add_argument("--modes", type=int, default=MODES)
    args = parser.parse_args()
    report = certify_negative_heads(
        args.matrix, expected_digest=args.sha256, cutoff=args.cutoff, half_width=args.half_width, modes=args.modes,
    )
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
