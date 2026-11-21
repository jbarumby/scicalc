"""Certify the raised-floor finite auxiliary and its threshold implications."""

import argparse
from fractions import Fraction
from hashlib import file_digest, sha256
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from NUMERIA_ES8P249.legendre_tail import modal_tail_bounds
from NUMERIA_ES8P251.gate_b.finite_head import _dyadic_fraction, load_aggregate
from NUMERIA_ES8P251.gate_b.finite_head_certify import assemble_finite_head, certify_center_weyl
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import fixed_beta
from NUMERIA_ES8P251.gate_c.spectral_diagnostic import _point, _rayleigh
from NUMERIA_ES8P253.prime7_tail import certify_four_prime_floor
from NUMERIA_ES8P257.p257_strip_taper import arctanh_lower, arctanh_upper

from .concentration import (
    CUTOFF, HALF_WIDTH, MODES, SECTORS, analytic_concentration, compare_concentration,
    gl60_error_budget, load_authenticated, matrix_rows,
)
from .independent_rayleigh import replay_rayleigh
from .strip_budget import LOSS, certify_strip_constants


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "NUMERIA_ES8P253/data"
OLD_AGGREGATE = DATA / "P253_CLOSE504.json.gz"
OLD_SHA = "20f86083413bb5fcf3a5faaa05aee035d0a9d015664ba0d5108013beeb21e914"
FLOOR = Fraction(1, 4)
LOWER_WINDOW = Fraction("0.9729550745276566525526763717215898648185")
DILATION_LIPSCHITZ = 21250


def _ball(value):
    return arb(value.numerator) / value.denominator


def _verify_tail():
    raw = (json.dumps(certify_four_prime_floor(), separators=(",", ":")) + "\n").encode()
    if raw != (DATA / "FOUR_PRIME_RAISED_FLOOR.json").read_bytes():
        raise ArithmeticError("P253 raised-floor tail replay differs from archived bytes")
    return sha256(raw).hexdigest()


def _old_matrices(path):
    with path.open("rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != OLD_SHA:
            raise ValueError("P253 weighted aggregate SHA-256 mismatch")
    matrices, _, digest = load_aggregate(path, half_width=HALF_WIDTH, cutoff=CUTOFF, modes=MODES)
    if digest != OLD_SHA:
        raise ValueError("P253 weighted aggregate changed while reading")
    return matrices


def update_head(head, gram, delta):
    """Apply delta times the high-pass projection compressed to the head."""
    if head.nrows() != gram.nrows() or head.ncols() != gram.ncols():
        raise ValueError("weighted head and concentration block dimensions differ")
    return head + delta * (arb_mat(head.nrows(), head.ncols(), [
        arb(i == j) for i in range(head.nrows()) for j in range(head.ncols())
    ]) - gram)


def _trial_after_failed_pivot(matrix, k):
    leading = arb_mat(k, k, [matrix[i, j].mid() for i in range(k) for j in range(k)])
    rhs = arb_mat(k, 1, [-matrix[i, k].mid() for i in range(k)])
    partial = [value.mid() for value in leading.solve(rhs).entries()] if k else []
    return partial + [arb(1)]


def _witness(matrix, k, errors):
    points = _trial_after_failed_pivot(matrix, k)
    vector = arb_mat(matrix.nrows(), 1, points + [arb(0)] * (matrix.nrows() - len(points)))
    quotient = _rayleigh(matrix, vector) + arb(0, _ball(sum(errors.values())).upper())
    if not quotient < 0:
        return None
    # Only rows on witness support are needed in the exact rational replay.
    leading = arb_mat(len(points), len(points), [
        matrix[i, j] for i in range(len(points)) for j in range(len(points))
    ])
    replay = replay_rayleigh(
        [_dyadic_fraction(point.man_exp()) for point in points],
        {"assembled": (matrix_rows(leading), 1)}, operator_errors=errors,
    )
    if replay["witness_strict_sign"] != "negative":
        raise ArithmeticError("the exact rational replay could not verify the negative witness")
    # Transfer by unitary dilation to a rational window below the threshold.
    if not LOWER_WINDOW < arctanh_lower(Fraction(3, 4), 180) / 2 < HALF_WIDTH:
        raise ArithmeticError("lower support window could not be verified")
    loss = DILATION_LIPSCHITZ * (HALF_WIDTH - LOWER_WINDOW)
    upper_lower_window = Fraction(replay["rayleigh_interval_exact"][1]) + loss
    return {
        "Arb_rayleigh_endpoints_dyadic": [_point(quotient.lower()), _point(quotient.upper())],
        "exact_rational_replay": replay,
        "coordinates_exact_dyadic": [_point(point) for point in points],
        "full_dimension": matrix.nrows(),
        "dilation_width_below_log7_over_2_exact": str(LOWER_WINDOW),
        "dilation_loss_exact": str(loss),
        "rayleigh_upper_below_threshold_exact": str(upper_lower_window),
        "negative_by_zero_extension_for_all_windows_beyond_threshold": upper_lower_window < 0,
    }


def _conditioning(matrix, error):
    """Bound kappa from below using one small-residual exact vector."""
    center = matrix.mid()
    n = center.nrows()
    vector = arb_mat(n, 1, [((i * i + 7 * i + 3) % 17) - 8 for i in range(n)])
    for _ in range(2):
        try:
            vector = center.solve(vector).mid()
        except ZeroDivisionError:
            return None
    residual = matrix * vector
    residual_sq = sum((value.abs_upper()**2 for value in residual.entries()), arb(0))
    norm_sq = sum((value**2 for value in vector.entries()), arb(0))
    small = ((residual_sq / norm_sq).sqrt() + _ball(error)).upper()
    large = (max(matrix[i, i].abs_lower() for i in range(n)) - _ball(error)).lower()
    if not (small > 0 and large > 0):
        return None
    return {
        "singular_value_min_upper_dyadic": _point(small),
        "operator_norm_lower_dyadic": _point(large),
        "condition_number_lower_dyadic": _point((large / small).lower()),
    }


def _audit_sector(old_quad, gram, sector, use_gl60):
    old, assembly = assemble_finite_head(old_quad, sector, half_width=HALF_WIDTH, cutoff=CUTOFF)
    delta = _ball(FLOOR) - fixed_beta()
    if not 0 < delta < _ball(FLOOR):
        raise ArithmeticError("unexpected floor ordering")
    delta_bound = _dyadic_fraction(delta.upper().man_exp())
    old_error = assembly["gl60_error_budget"]["row_sum_bound"]
    gram_error = gl60_error_budget(sector)["row_sum_bound"] if use_gl60 else Fraction(0)
    errors = {"weighted_GL60": old_error, "concentration_GL60_scaled": delta_bound * gram_error}
    error = sum(errors.values())
    matrix = update_head(old, gram, delta)
    positivity = certify_center_weyl(matrix, error)
    result = {
        "validation": positivity,
        "operator_error_budget_exact": {name: str(value) for name, value in errors.items()},
        "condition": _conditioning(matrix, error),
    }
    if not positivity["certified_positive"]:
        k = positivity.get("first_inconclusive_pivot_index")
        result["negative_witness"] = None if k is None else _witness(matrix, k, errors)
    return result


def _schur(sectors):
    if not all(result["validation"]["certified_positive"] for result in sectors.values()):
        return None
    results = {}
    for parity, sector in enumerate(SECTORS):
        coupling_squared, deviation = modal_tail_bounds(CUTOFF, MODES, parity, half_width=HALF_WIDTH)
        head = _dyadic_fraction(sectors[sector]["validation"]["margin_lower_dyadic"])
        complement = FLOOR - deviation
        margin = min(head, complement) / 2
        if not margin > 0 or (head - margin) * (complement - margin) <= coupling_squared:
            return None
        results[sector] = str(margin)
    margin = min(map(Fraction, results.values()))
    room = 2 * HALF_WIDTH - arctanh_upper(Fraction(3, 4), 180)
    if room <= 0:
        raise ArithmeticError("no verified crossing room")
    strip = min(room / 2, margin / (2 * LOSS))
    return {
        "auxiliary_sector_margin_lower_exact": results,
        "auxiliary_margin_lower_exact": str(margin),
        "Weil_half_width_exact_expression": "log(7)/2 + delta/2",
        "delta_exact": str(strip),
        "Weil_margin_lower_exact": str(margin - LOSS * strip),
    }


def certify(old_path, gram_path=None, gram_sha=None):
    if (gram_path is None) != (gram_sha is None):
        raise ValueError("provide both the GL60 aggregate and its SHA-256")
    old = _old_matrices(Path(old_path))
    tail_sha = _verify_tail()
    endpoint = analytic_concentration()
    checked = None
    grams = endpoint
    if gram_path is not None:
        grams = load_authenticated(gram_path, gram_sha)
        checked = compare_concentration(grams, endpoint)
    if not arctanh_upper(Fraction(3, 4), 180) / 2 < HALF_WIDTH:
        raise ArithmeticError("tested width does not lie above the first new prime threshold")
    with ctx.workprec(1024):
        sectors = {s: _audit_sector(q, g, s, checked is not None)
                   for s, q, g in zip(SECTORS, old, grams, strict=True)}
    return {
        "schema": "NUMERIA.ES8P258.RAISED_AUXILIARY.1",
        "scope": "four-prime truncated auxiliary at cutoff 504",
        "L_exact": str(HALF_WIDTH),
        "cutoff": CUTOFF,
        "modes_per_sector": MODES,
        "beta_raised_exact": str(FLOOR),
        "P253_weighted_GL60_sha256": OLD_SHA,
        "P253_raised_tail_replay_sha256": tail_sha,
        "G504_GL60_sha256": gram_sha,
        "G504_GL60_endpoint_crosscheck": checked,
        "sectors": sectors,
        "prime7_strip_conditional_constants": certify_strip_constants(),
        "model_and_local_Weil_margins": _schur(sectors),
        "cofinal_extension_established": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, default=OLD_AGGREGATE)
    parser.add_argument("--gram", type=Path)
    parser.add_argument("--gram-sha")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(certify(args.old, args.gram, args.gram_sha), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
