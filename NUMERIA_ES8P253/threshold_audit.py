"""Replay the four-prime model and certify a window across the log(7)/2 threshold."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from NUMERIA_ES8P249.legendre_tail import ARCHIVED_CUTOFF, four_prime_model_half_width
from NUMERIA_ES8P250.adaptive_signed_tail import run as replay_four_prime_tail
from NUMERIA_ES8P251.gate_b.finite_head import (
    _dyadic_fraction,
    _schur_transfer,
    aggregate_chunks,
    load_aggregate,
    write_aggregate,
)
from NUMERIA_ES8P251.gate_b.finite_head_certify import certify_finite_head
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import CUTOFF, MODES, SECTORS, fixed_beta
from NUMERIA_ES8P251.p251_transfer_gate import PINNED_P250_SHA256


AUX_WIDTH = Fraction("0.972955074527656652552676371721590")
OMEGA = 4 * ARCHIVED_CUTOFF
LOSS_PER_STRIP_WIDTH = 2520


def crossing_bounds(aux_floor, *, half_width=AUX_WIDTH):
    """Verify the constants in the boundary-strip form bound."""
    aux_floor = Fraction(aux_floor)
    if aux_floor <= 0:
        raise ValueError("the auxiliary four-prime model needs positive coercivity")
    half_width = four_prime_model_half_width(half_width)
    with ctx.workprec(256):
        shift = arb(7).log()
        kappa = shift / arb(7).sqrt()
        aux_span = arb(2 * half_width.numerator) / half_width.denominator
        headroom = _dyadic_fraction((aux_span - shift).lower().man_exp())
        if headroom <= 0:
            raise ValueError("the auxiliary width must reach strictly beyond log(7)/2")
        strip_width = min(aux_floor / (2 * LOSS_PER_STRIP_WIDTH), headroom / 2)
        span = shift + arb(strip_width.numerator) / strip_width.denominator
        if not (span < aux_span < arb(8).log() and aux_span < 2 * shift):
            raise ArithmeticError("the new shift is outside the two-strip regime")
        reserve = arb(4).log() + arb(1) / ARCHIVED_CUTOFF - arb(1) / OMEGA
        if not (0 < fixed_beta() < kappa < arb(3) / 4 and arb.pi() > 3 and reserve > arb(4) / 3):
            raise ArithmeticError("the interval proof of the strip constants failed")
    crossing_floor = aux_floor - LOSS_PER_STRIP_WIDTH * strip_width
    return {
        "four_prime_auxiliary_floor_exact": str(aux_floor),
        "auxiliary_span_headroom_lower_exact": str(headroom),
        "cutoff_for_absorption": OMEGA,
        "prime7_amplitude_strict_upper_exact": "3/4",
        "four_prime_high_tail_floor_strict_lower_exact": "4/3",
        "strip_width_exact": str(strip_width),
        "half_width_increment_above_log7_over_2_exact": str(strip_width / 2),
        "loss_per_strip_width_upper": LOSS_PER_STRIP_WIDTH,
        "Weil_floor_lower_exact": str(crossing_floor),
        "additional_high_frequency_floor_lower_exact": "1/12",
    }


def certify_crossing(matrix, *, expected_digest, cutoff=350, half_width=AUX_WIDTH, modes=MODES):
    if not CUTOFF <= cutoff <= ARCHIVED_CUTOFF:
        raise ValueError("the retained band must end between the P250 and archived cutoffs")
    matrices, chunks, digest = load_aggregate(matrix, cutoff=cutoff, half_width=half_width, modes=modes)
    if digest != expected_digest:
        raise ArithmeticError("P253 auxiliary GL60 aggregate SHA-256 mismatch")
    tail = replay_four_prime_tail()
    tail_bytes = (json.dumps(tail, indent=2, ensure_ascii=False) + "\n").encode()
    tail_digest = sha256(tail_bytes).hexdigest()
    if (tail_digest != PINNED_P250_SHA256
            or tail["certified_positive_frequency_range"] != [CUTOFF, "infinity"]
            or tail["beta_frozen_at"] != ARCHIVED_CUTOFF):
        raise ArithmeticError("four-prime tail certificate mismatch")

    sectors = {
        sector: certify_finite_head(quadrature, sector, cutoff=cutoff, half_width=half_width)
        for sector, quadrature in zip(SECTORS, matrices, strict=True)
    }
    schur = _schur_transfer(sectors, cutoff=cutoff, half_width=half_width, modes=modes)
    crossing = None
    if schur["certified_for_reduced_model"]:
        crossing = crossing_bounds(schur["model_floor_lower_exact"], half_width=half_width)
    return {
        "schema": "NUMERIA.ES8P253.PRIME7_THRESHOLD.1",
        "auxiliary_half_width_exact": str(half_width),
        "auxiliary_prime_powers": [2, 3, 4, 5],
        "prime_powers_at_new_endpoint": [2, 3, 4, 5, 7],
        "cutoff": cutoff,
        "modes_per_sector": modes,
        "input_GL60_aggregate_sha256": digest,
        "input_GL60_chunks": chunks,
        "signed_tail_P250_recomputed_sha256": tail_digest,
        "auxiliary_sectors": sectors,
        "auxiliary_modal_Schur": schur,
        "threshold_crossing": crossing,
        "Weil_test_function_class": "complex C_c^infty(R) supported in [-L,L]",
        "cofinal_extension_established": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunks", nargs="+", type=Path)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--half-width", type=Fraction, default=AUX_WIDTH)
    parser.add_argument("--cutoff", type=int, default=350)
    parser.add_argument("--modes", type=int, default=MODES)
    parser.add_argument("--sha256", required=True, help="expected compressed GL60 aggregate digest")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.chunks:
        matrices, records = aggregate_chunks(
            args.chunks, cutoff=args.cutoff, half_width=args.half_width, modes=args.modes,
        )
        write_aggregate(args.matrix, matrices, records, cutoff=args.cutoff, half_width=args.half_width, modes=args.modes)
    result = certify_crossing(
        args.matrix, expected_digest=args.sha256, cutoff=args.cutoff, half_width=args.half_width, modes=args.modes,
    )
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
