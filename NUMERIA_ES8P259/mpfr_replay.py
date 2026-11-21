"""Independently verify the P259 cut-590 finite heads with directed MPFR."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from flint import ctx

from NUMERIA_ES8P251.gate_b.finite_head import _dyadic_fraction
from NUMERIA_ES8P251.gate_b.mpfr_audit import (
    P242_SOURCE, P242_SOURCE_SHA256, _export_centers, _read_trace, _weyl_margin,
)
from NUMERIA_ES8P257.p257_strip_taper import arctanh_lower, arctanh_upper
from NUMERIA_ES8P258.audit import CUTOFF, HALF_WIDTH, MODES, OLD_AGGREGATE, OLD_SHA, SECTORS, _old_matrices, _verify_tail
from NUMERIA_ES8P258.concentration import analytic_concentration

from .banded_weil import LOSS, _round_up_modal_bound, assemble_banded_head, certify_analytic_constants, chebyshev_tail_bounds


HIGH_CUTOFF = 590
PRECISION = 1024


def replay(compiler="cc"):
    source_bytes = P242_SOURCE.read_bytes()
    if sha256(source_bytes).hexdigest() != P242_SOURCE_SHA256:
        raise ValueError("P242 C verifier SHA-256 mismatch")
    constants = certify_analytic_constants(HIGH_CUTOFF)
    tail_sha = _verify_tail()
    old = _old_matrices(OLD_AGGREGATE)
    g504 = analytic_concentration(cutoff=CUTOFF)
    g590 = analytic_concentration(cutoff=HIGH_CUTOFF, bits=2048)

    sectors = {}
    with TemporaryDirectory(prefix="p259-mpfr-") as directory:
        path = Path(directory)
        source, binary = path / "p242.c", path / "interval-ldl"
        source.write_bytes(source_bytes)
        subprocess.run([compiler, "-O3", "-Wall", str(source), "-o", str(binary), "-l:libmpfr.so.6", "-lgmp"], check=True)
        for parity, (sector, weighted, gram, extended) in enumerate(zip(SECTORS, old, g504, g590, strict=True)):
            centers, trace = path / f"{sector}.centers", path / f"{sector}.trace"
            with ctx.workprec(PRECISION):
                _, matrix, analytic = assemble_banded_head(weighted, gram, extended, sector)
                radius = _export_centers(matrix, centers)
            subprocess.run([str(binary), str(centers), str(trace), str(PRECISION), str(MODES)], check=True)
            output = trace.read_bytes()
            pivot, inverse_squared = _read_trace(output, MODES)
            head = _weyl_margin(_dyadic_fraction(inverse_squared), radius + analytic)
            cross, deviation = chebyshev_tail_bounds(HIGH_CUTOFF, MODES, parity)
            complement = Fraction(1, 2) - _round_up_modal_bound(deviation)
            margin = min(head, complement) / 2
            if not (margin > 0 and (head - margin) * (complement - margin) > _round_up_modal_bound(cross)):
                raise ArithmeticError(f"independent Schur estimate failed for {sector}")
            sectors[sector] = {
                "positive_pivots": MODES,
                "smallest_pivot_lower_dyadic": pivot,
                "inverse_frobenius_squared_upper_dyadic": inverse_squared,
                "centers_sha256": sha256(centers.read_bytes()).hexdigest(),
                "trace_sha256": sha256(output).hexdigest(),
                "finite_head_margin_exact": str(head),
                "auxiliary_margin_exact": str(margin),
            }

    lower = arctanh_lower(Fraction(3, 4), 180)
    upper = arctanh_upper(Fraction(3, 4), 180)
    margin = min(Fraction(record["auxiliary_margin_exact"]) for record in sectors.values())
    adjusted = margin - LOSS * (2 * HALF_WIDTH - lower)
    if not (upper < 2 * HALF_WIDTH < 3 * arctanh_lower(Fraction(1, 3), 4) and adjusted > 0):
        raise ArithmeticError("independent edge budget did not certify the target window")

    return {
        "schema": "NUMERIA.ES8P259.DIRECTED_MPFR_REPLAY.1",
        "precision_bits": PRECISION,
        "old_weighted_GL60_sha256": OLD_SHA,
        "old_four_prime_tail_replay_sha256": tail_sha,
        "verifier_sha256": P242_SOURCE_SHA256,
        "high_band_start": HIGH_CUTOFF,
        "modes_per_sector": MODES,
        "sectors": sectors,
        "strip_loss_integer_upper": constants["strip_loss_integer_upper"],
        "Weil_half_width_exact": str(HALF_WIDTH),
        "Weil_margin_lower_exact": str(adjusted),
        "cofinal_extension_established": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", default="cc")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(replay(args.cc), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
