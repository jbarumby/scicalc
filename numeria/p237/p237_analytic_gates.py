#!/usr/bin/env python3
"""Exact rational geometry gates for the P237 rho=3 Bernstein ellipses.

This module derives the E3 panel bounds for the digamma and spherical-Bessel
arguments using Fraction arithmetic. Special-function identities and the
Bernstein/Gauss theorems are classical analytic inputs, not formalized here.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def make_gates() -> dict:
    L = Q(161, 200)
    rho = Q(3)
    panel_scale = Q(1, 4)
    last_panel = 1007

    # Bernstein ellipse E_rho semiaxes, then the t=(center+s/4) map.
    ellipse_a = (rho + 1 / rho) / 2
    ellipse_b = (rho - 1 / rho) / 2
    t_real_semiaxis = panel_scale * ellipse_a
    t_imag_semiaxis = panel_scale * ellipse_b
    panel_center = Q(2 * last_panel + 1, 4)
    t_real_max = panel_center + t_real_semiaxis
    t_imag_max = t_imag_semiaxis

    # Both holomorphic digamma arguments z_±=1/4±it/2 stay in Re(z)>0.
    digamma_real_min = Q(1, 4) - t_imag_max / 2
    digamma_real_max = Q(1, 4) + t_imag_max / 2
    digamma_imag_max = t_real_max / 2
    digamma_mod_l1_upper = digamma_real_max + digamma_imag_max

    # Rectangular enclosure of Lt. Its squared Euclidean radius is exact.
    bessel_real_max = L * t_real_max
    bessel_imag_max = L * t_imag_max
    bessel_mod_sq_upper = bessel_real_max**2 + bessel_imag_max**2

    checks = {
        "ellipse_real_semiaxis_5_over_3": ellipse_a == Q(5, 3),
        "ellipse_imag_semiaxis_4_over_3": ellipse_b == Q(4, 3),
        "t_real_max_3025_over_6": t_real_max == Q(3025, 6),
        "t_imag_max_1_over_3": t_imag_max == Q(1, 3),
        "digamma_real_min_1_over_12": digamma_real_min == Q(1, 12),
        "digamma_real_max_5_over_12": digamma_real_max == Q(5, 12),
        "digamma_modulus_lt_253": digamma_mod_l1_upper < 253,
        "bessel_real_bound_19481_over_48": bessel_real_max == Q(19481, 48),
        "bessel_imag_bound_161_over_600": bessel_imag_max == Q(161, 600),
        "bessel_imag_lt_one_half": bessel_imag_max < Q(1, 2),
        "bessel_modulus_lt_406": bessel_mod_sq_upper < 406**2,
        "uniform_exponential_bound_lt_2": 1 / (1 - bessel_imag_max) < 2,
        "panel_count_1008": 2 * 504 == 1008,
        "even_mode_count_377": len(range(0, 753, 2)) == 377,
        "odd_mode_count_377": len(range(1, 754, 2)) == 377,
    }
    need(all(checks.values()), "an exact rational E3 geometry gate failed")

    return {
        "schema": "NUMERIA.ES8P237.ANALYTIC_GATES.2",
        "status": "exact_rational_geometry_verified_classical_analytic_identities_imported",
        "L": str(L),
        "Tsharp": 504,
        "Bernstein_rho": 3,
        "panels": 1008,
        "t_real_abs_upper": str(t_real_max),
        "t_imag_abs_upper": str(t_imag_max),
        "digamma_real_range": [str(digamma_real_min), str(digamma_real_max)],
        "digamma_imag_abs_upper": str(digamma_imag_max),
        "digamma_modulus_l1_upper": str(digamma_mod_l1_upper),
        "bessel_Lt_real_abs_upper": str(bessel_real_max),
        "bessel_Lt_imag_abs_upper": str(bessel_imag_max),
        "bessel_Lt_modulus_squared_upper": str(bessel_mod_sq_upper),
        "bessel_Lt_modulus_upper": 406,
        "uniform_exp_upper_from_geometric_series": str(1 / (1 - bessel_imag_max)),
        "checks": checks,
        "analytic_inputs": [
            "digamma partial-fraction identity and conjugation",
            "spherical-Bessel plane-wave and Rodrigues integral identities",
            "Legendre orthogonality",
            "Chebyshev coefficient bound on Bernstein ellipses",
            "Gauss-Legendre exactness and positive weights",
        ],
        "Weil_Zhu_identity_verified": False,
        "infinite_dimensional_Weil_positivity_proved": False,
        "RH_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = make_gates()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"PASS exact E3 geometry gates: {len(result['checks'])}")


if __name__ == "__main__":
    main()
