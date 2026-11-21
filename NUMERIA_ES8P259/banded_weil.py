"""Test the Weil auxiliary with piecewise constant floors 1/4 and 1/2."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P249.interval_audit import _arb_prime_terms
from NUMERIA_ES8P249.legendre_tail import modal_tail_bounds
from NUMERIA_ES8P249.prime_operator import P236_PRIME_POWERS
from NUMERIA_ES8P251.gate_b.finite_head import _dyadic_fraction
from NUMERIA_ES8P251.gate_b.finite_head_certify import assemble_finite_head, certify_center_weyl
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import fixed_beta
from NUMERIA_ES8P251.gate_c.spectral_diagnostic import _point, _rayleigh
from NUMERIA_ES8P253.prime7_tail import _positive_cover
from NUMERIA_ES8P257.p257_strip_taper import (
    HIGH_REMAINDER, KAPPA_UP, PI_LO_TARGET, arctanh_lower, arctanh_upper, prime_7_constants,
)
from NUMERIA_ES8P258.audit import (
    CUTOFF, FLOOR, HALF_WIDTH, MODES, OLD_AGGREGATE, OLD_SHA, SECTORS, _ball,
    _old_matrices, _trial_after_failed_pivot, _verify_tail, update_head,
)
from NUMERIA_ES8P258.concentration import analytic_concentration, matrix_rows
from NUMERIA_ES8P258.independent_rayleigh import replay_rayleigh


HIGH_CUTOFF = 850
STEP = Fraction(1, 4)
OMEGA = 3360
LOSS = 1839
ROOT = Path(__file__).resolve().parents[1]


def certify_analytic_constants(high_cutoff=HIGH_CUTOFF):
    prime_7_constants()
    log850_over504 = Fraction(2 * (HIGH_CUTOFF - CUTOFF), HIGH_CUTOFF + CUTOFF)
    if not log850_over504 > Fraction(1, 2):
        raise ArithmeticError("insufficient Binet bound above 850")
    adaptive_certificate = None
    if high_cutoff == 590:
        with ctx.workprec(192):
            adaptive_certificate = _positive_cover(
                _arb_prime_terms(P236_PRIME_POWERS), arb(1) / 2, high_cutoff, HIGH_CUTOFF, arb.pi(),
            )
    elif high_cutoff != HIGH_CUTOFF:
        raise ValueError("the supported high-band cutoffs are 590 and 850")
    tau = Fraction(13, 10)
    reserve = arctanh_lower(Fraction(17, 23), 8) - Fraction(1, 2)
    loss = (OMEGA + 1) * KAPPA_UP * tau / (PI_LO_TARGET * (tau - KAPPA_UP))
    if not (reserve > tau + HIGH_REMAINDER and loss < LOSS):
        raise ArithmeticError("insufficient boundary-strip reserve")
    return {
        "high_band_start": high_cutoff,
        "high_band_adaptive_certificate": adaptive_certificate,
        "high_band_floor_exact": "1/2",
        "high_band_log_ratio_strict_lower_exact": str(log850_over504),
        "absorption_cutoff": OMEGA,
        "higher_tail_reserve_strict_lower_exact": str(reserve),
        "strip_loss_rational_upper_exact": str(loss),
        "strip_loss_integer_upper": LOSS,
        "high_frequency_margin_exact": str(HIGH_REMAINDER),
    }


def effective_head_step(step, reserved):
    """Price of using part of the complement to control the positive update."""
    if step <= 0 or reserved <= 0:
        raise ValueError("both the step and the reserved complement must be positive")
    return step * reserved / (step + reserved)


def _negative_witness(matrix, check, operator_error):
    k = check.get("first_inconclusive_pivot_index")
    if k is None:
        return None
    trial = _trial_after_failed_pivot(matrix, k)
    leading = arb_mat(len(trial), len(trial), [
        matrix[i, j] for i in range(len(trial)) for j in range(len(trial))
    ])
    report = replay_rayleigh(
        [_dyadic_fraction(v.man_exp()) for v in trial],
        {"raised_finite_head": (matrix_rows(leading), 1)},
        operator_errors={"base_weighted_GL60": operator_error},
    )
    return report if report["witness_strict_sign"] == "negative" else None


def _past_witness(matrix, sector, operator_error):
    if sector != "odd":
        return None
    summary = json.loads((ROOT / "NUMERIA_ES8P258/P258_RESULTADO.json").read_text(encoding="utf-8"))
    points = [arb((int(m), int(e))) for m, e in summary["odd_witness_nonpadding_coordinates_exact_dyadic"]]
    vector = arb_mat(MODES, 1, points + [arb(0)] * (MODES - len(points)))
    q = _rayleigh(matrix, vector) + arb(0, _ball(operator_error).upper())
    return [_point(q.lower()), _point(q.upper())]


def _round_up_modal_bound(value):
    # Keep the proof rational and the serialized denominator bounded.
    scale = 10**110
    return Fraction(-(-value.numerator * scale // value.denominator), scale)


def _sector(weighted, gram504, gram_high, sector, parity):
    original, metadata = assemble_finite_head(weighted, sector, half_width=HALF_WIDTH, cutoff=CUTOFF)
    error = metadata["gl60_error_budget"]["row_sum_bound"]
    base = update_head(original, gram504, _ball(FLOOR) - fixed_beta())
    raised = update_head(base, gram_high, _ball(STEP))
    check = certify_center_weyl(raised, error)
    result = {
        "finite_head": check,
        "old_P258_witness_rayleigh_in_new_head_dyadic": _past_witness(raised, sector, error),
        "continuous_model_margin_exact": None,
        "block_transfers": [],
    }
    if not check["certified_positive"]:
        result["finite_negative_witness"] = _negative_witness(raised, check, error)
        return result

    coupling_sq, deviation = modal_tail_bounds(CUTOFF, MODES, parity, half_width=HALF_WIDTH)
    coupling_sq = _round_up_modal_bound(coupling_sq)
    tail_lower = FLOOR - _round_up_modal_bound(deviation)
    for reserved in (Fraction(1, 8), Fraction(6, 25), Fraction(249, 1000), Fraction(24999, 100000)):
        complement_lower = tail_lower - reserved
        if not complement_lower > 0:
            raise ArithmeticError("no remainder of complement margin")
        factor = effective_head_step(STEP, reserved)
        conservative = update_head(base, gram_high, _ball(factor))
        test = certify_center_weyl(conservative, error)
        record = {
            "reserved_complement_exact": str(reserved),
            "effective_update_exact": str(factor),
            "head_test": test,
            "complement_lower_exact": str(complement_lower),
            "coupling_squared_upper_exact": str(coupling_sq),
        }
        result["block_transfers"].append(record)
        if test["certified_positive"]:
            head = _dyadic_fraction(test["margin_lower_dyadic"])
            margin = min(head, complement_lower) / 2
            if (head - margin) * (complement_lower - margin) > coupling_sq:
                result["continuous_model_margin_exact"] = str(margin)
                break
    return result


def certify(high_cutoff=HIGH_CUTOFF):
    constants = certify_analytic_constants(high_cutoff)
    tail_sha = _verify_tail()
    weighted = _old_matrices(OLD_AGGREGATE)
    gram504 = analytic_concentration(cutoff=CUTOFF)
    gram_high = analytic_concentration(cutoff=high_cutoff, bits=2048)
    with ctx.workprec(1024):
        sectors = {s: _sector(q, a, b, s, parity)
                   for parity, (s, q, a, b) in enumerate(zip(SECTORS, weighted, gram504, gram_high, strict=True))}
    crossing = None
    if all(record["continuous_model_margin_exact"] is not None for record in sectors.values()):
        margin = min(Fraction(record["continuous_model_margin_exact"]) for record in sectors.values())
        room = 2 * HALF_WIDTH - arctanh_upper(Fraction(3, 4), 180)
        if not room > 0:
            raise ArithmeticError("no certified crossing room")
        strip = min(room / 2, margin / (2 * LOSS))
        crossing = {
            "delta_exact": str(strip),
            "new_Weil_half_width_expression": "log(7)/2 + delta/2",
            "Weil_strict_margin_lower_exact": str(margin - LOSS * strip),
        }
    return {
        "schema": "NUMERIA.ES8P259.WEIL_BANDED.1",
        "four_prime_weighted_GL60_input_sha256": OLD_SHA,
        "floor_1_4_tail_replay_sha256": tail_sha,
        "auxiliary_L_exact": str(HALF_WIDTH),
        "modes_per_sector": MODES,
        "band_and_strip": constants,
        "sectors": sectors,
        "Weil_crossing": crossing,
        "cofinal_extension_established": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--high-cutoff", type=int, default=HIGH_CUTOFF, choices=(590, HIGH_CUTOFF))
    args = parser.parse_args()
    args.out.write_text(json.dumps(certify(args.high_cutoff), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
