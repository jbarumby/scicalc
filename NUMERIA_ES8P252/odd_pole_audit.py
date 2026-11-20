"""Certify the P252 finite odd head by a scalar test for its negative pole."""

import argparse
from fractions import Fraction
import json
from operator import index
from pathlib import Path

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P251.gate_b.finite_head import _point, load_aggregate
from NUMERIA_ES8P251.gate_b.finite_head_certify import (
    _max_radius_row_sum,
    certify_center_weyl,
    gl60_error_budget,
)
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import CUTOFF, MODES, fixed_beta

from .window_audit import HALF_WIDTH


P252_GL60_SHA256 = "85cccdafd673222792b59ed1e1c0d25f3c39586fbac808958c4af9619c65bbe7"
REFERENCE_MARGIN = Fraction("1e-23")
SCHEMA = "NUMERIA.ES8P252.ODD_POLE_SCHUR.1"


def certify_odd_pole(path, *, margin=REFERENCE_MARGIN, precision_bits=1024):
    """Check A−2ppᵀ ≻ mI using a lower bound M ≤ A−mI and a scalar solve.

    A is the pole-free continuous Legendre head for L=19/20. The interval
    quadrature perturbation and the GL60 analytical budget are subtracted
    from the pole-free midpoint before certifying M. For M positive definite,
    M−2ppᵀ ≻ 0 is equivalent to 1−2pᵀ M⁻¹p > 0 by congruence with M⁻¹ᐟ².
    """
    margin, precision_bits = Fraction(margin), index(precision_bits)
    if margin <= 0 or precision_bits < 128:
        raise ValueError("expected positive margin and at least 128 bits")
    (_, odd), _, digest = load_aggregate(path, half_width=HALF_WIDTH)
    if digest != P252_GL60_SHA256:
        raise ValueError("P252 quadrature aggregate SHA-256 mismatch")

    with ctx.workprec(precision_bits):
        half_width = arb(HALF_WIDTH.numerator) / HALF_WIDTH.denominator
        requested_upper = (arb(margin.numerator) / margin.denominator).upper()
        beta_lower = fixed_beta().lower()
        radius_upper = _max_radius_row_sum(odd)
        gl_error = gl60_error_budget("odd", MODES, CUTOFF, half_width=HALF_WIDTH)["row_sum_bound"]
        gl_upper = (arb(gl_error.numerator) / gl_error.denominator).upper()

        guarded_base = odd.mid()
        for k in range(MODES):
            guarded_base[k, k] += beta_lower - requested_upper - radius_upper - gl_upper
        base_proof = certify_center_weyl(guarded_base)
        if not base_proof["certified_positive"]:
            raise ArithmeticError("the shifted pole-free odd head is not certified positive")

        z = half_width / 2
        prefactor = (arb.pi() / (2 * z)).sqrt()
        pole = arb_mat(MODES, 1, [
            (2 * half_width * (2 * degree + 1)).sqrt()
            * prefactor * z.bessel_i(arb(2 * degree + 1) / 2)
            for degree in range(1, 2 * MODES, 2)
        ])
        scalar_guard = 1 - 2 * (pole.transpose() * guarded_base.solve(pole))[0, 0]
        if not scalar_guard.is_finite() or not scalar_guard.lower() > 0:
            raise ArithmeticError("the interval scalar test does not certify the negative pole")

        return {
            "schema": SCHEMA,
            "input_GL60_aggregate_sha256": digest,
            "L": str(HALF_WIDTH),
            "cutoff": CUTOFF,
            "sector": "odd",
            "modes": MODES,
            "precision_bits": precision_bits,
            "finite_head_margin_strict_lower_exact": str(margin),
            "quadrature_radius_row_upper_dyadic": _point(radius_upper),
            "gl60_analytic_row_upper_exact": str(gl_error),
            "shifted_pole_free_positive_pivots": base_proof["positive_pivot_count"],
            "shifted_pole_free_margin_lower_dyadic": base_proof["margin_lower_dyadic"],
            "one_minus_twice_pole_inverse_pole_endpoints_dyadic": [
                _point(scalar_guard.lower()), _point(scalar_guard.upper())
            ],
            "finite_odd_head_certified": True,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="P252 aggregate GL60 matrix")
    parser.add_argument("--out", type=Path, required=True, help="odd-pole audit JSON")
    parser.add_argument("--margin", type=Fraction, default=REFERENCE_MARGIN)
    parser.add_argument("--precision-bits", type=int, default=1024)
    args = parser.parse_args()
    report = certify_odd_pole(args.matrix, margin=args.margin, precision_bits=args.precision_bits)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
