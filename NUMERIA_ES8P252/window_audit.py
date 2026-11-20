"""Certify a wider Legendre window with the same signed prime symbol."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from NUMERIA_ES8P249.legendre_tail import fixed_symbol_half_width, modal_tail_bounds
from NUMERIA_ES8P250.adaptive_signed_tail import run as replay_signed_tail
from NUMERIA_ES8P251.gate_b.finite_head import _schur_transfer, aggregate_chunks, load_aggregate, write_aggregate
from NUMERIA_ES8P251.gate_b.finite_head_certify import certify_finite_head
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import CUTOFF, MODES, SECTORS
from NUMERIA_ES8P251.p251_transfer_gate import PINNED_P250_SHA256


HALF_WIDTH = Fraction(19, 20)


def certify_window(matrix, *, half_width=HALF_WIDTH, modes=MODES, precision_bits=1024):
    half_width = fixed_symbol_half_width(half_width)
    quadratures, chunks, digest = load_aggregate(matrix, modes=modes, half_width=half_width)
    sectors = {
        name: certify_finite_head(quadrature, name, precision_bits=precision_bits, half_width=half_width)
        for name, quadrature in zip(SECTORS, quadratures, strict=True)
    }
    tail = replay_signed_tail()
    tail_wire = (json.dumps(tail, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    if (sha256(tail_wire).hexdigest() != PINNED_P250_SHA256
            or tail["certified_positive_frequency_range"] != [CUTOFF, "infinity"]
            or tail["beta_frozen_at"] != 504):
        raise ArithmeticError("P250 fixed-beta signed tail failed to replay")

    schur = _schur_transfer(sectors, modes=modes, half_width=half_width)
    modal = {}
    for parity, name in enumerate(SECTORS):
        coupling_squared, deviation = modal_tail_bounds(CUTOFF, modes, parity, half_width=half_width)
        modal[name] = {
            "coupling_norm_strictly_below_1e_minus_54": coupling_squared < Fraction(1, 10**108),
            "tail_deviation_norm_strictly_below_1e_minus_117": deviation < Fraction(1, 10**117),
        }
    return {
        "schema": "NUMERIA.ES8P252.WIDER_WINDOW.1",
        "L": str(half_width),
        "cutoff": CUTOFF,
        "beta_reference_cutoff": 504,
        "panels_per_sector": 2 * CUTOFF,
        "GL60_nodes_per_sector": 120 * CUTOFF,
        "modes_per_sector": modes,
        "input_GL60_aggregate_sha256": digest,
        "input_GL60_chunks": chunks,
        "prime_powers": [2, 3, 4, 5],
        "signed_tail_P250_recomputed_sha256": PINNED_P250_SHA256,
        "sectors": sectors,
        "both_finite_heads_positive": all(sector["certified_positive"] for sector in sectors.values()),
        "modal_tail_budgets": modal,
        "schur_transfer": schur,
        "Weil_test_function_class": "complex C_c^infty(R) supported in [-L,L]",
        "reuses_Guinand_Weil_identification_from_P251_Gate_C": True,
        "cofinal_extension_established": False,
        "Riemann_hypothesis_proved": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunks", nargs="+", type=Path)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--half-width", type=Fraction, default=HALF_WIDTH)
    parser.add_argument("--modes", type=int, default=MODES)
    parser.add_argument("--precision-bits", type=int, default=1024)
    args = parser.parse_args()
    if args.chunks is not None:
        matrices, records = aggregate_chunks(args.chunks, modes=args.modes, half_width=args.half_width)
        write_aggregate(args.matrix, matrices, records, modes=args.modes, half_width=args.half_width)
    report = certify_window(
        args.matrix, modes=args.modes, half_width=args.half_width, precision_bits=args.precision_bits,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for name, certificate in report["sectors"].items():
        print(name, certificate["status"], flush=True)


if __name__ == "__main__":
    main()
