"""Check the P258 result archive and replay the even head in C/MPFR."""

import argparse
from fractions import Fraction
from hashlib import file_digest, sha256
import json
from math import isqrt
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from flint import acb, arb, ctx

from NUMERIA_ES8P251.gate_b.finite_head import _dyadic_fraction
from NUMERIA_ES8P251.gate_b.finite_head_certify import assemble_finite_head, certify_center_weyl
from NUMERIA_ES8P251.gate_b.finite_head_quadrature import _prime_terms, fixed_beta
from NUMERIA_ES8P251.gate_b.mpfr_audit import P242_SOURCE, P242_SOURCE_SHA256, _export_centers, _read_trace

from .audit import CUTOFF, FLOOR, HALF_WIDTH, MODES, OLD_AGGREGATE, _ball, _old_matrices, _witness, update_head
from .concentration import gl60_error_budget, load_authenticated


PINS = {
    "P258_ENDPOINT.json": "417e7447ab84af43fcb3af07444ac58a09e0d6ed9b4e006a719b7fda701930c1",
    "P258_G504.json.gz": "5826325eefdc4db73dcc1676c1dd8cdee23dd60d263fa32afb454e2202e6b041",
    "P258_GL60.json": "1aa765357cfe1fa46f864d5b57fc407157375e75b3fa8b23b77f0db5c2e4985c",
}


def _load_inputs(replay):
    for name, expected in PINS.items():
        with (replay / name).open("rb") as source:
            if file_digest(source, "sha256").hexdigest() != expected:
                raise ValueError(f"incorrect SHA-256 for {name}")
    result = json.loads((replay / "P258_GL60.json").read_text(encoding="utf-8"))
    gram = load_authenticated(replay / "P258_G504.json.gz", PINS["P258_G504.json.gz"])
    return result, gram


def _assemble(quadrature, gram, sector, floor):
    old, budget = assemble_finite_head(quadrature, sector, cutoff=CUTOFF, half_width=HALF_WIDTH)
    delta = _ball(floor) - fixed_beta()
    e = _dyadic_fraction(delta.abs_upper().man_exp())
    return update_head(old, gram, delta), {
        "weighted_GL60": budget["gl60_error_budget"]["row_sum_bound"],
        "concentration_GL60_scaled": e * gl60_error_budget(sector)["row_sum_bound"],
    }


def _mpfr_even(quadrature, gram):
    with ctx.workprec(1024):
        even, errors = _assemble(quadrature, gram, "even", FLOOR)
        with TemporaryDirectory(prefix="p258-mpfr-") as folder:
            work = Path(folder)
            source, binary, centers, trace = (work / f for f in ("ldl.c", "ldl", "centers.txt", "trace.txt"))
            source_bytes = P242_SOURCE.read_bytes()
            if sha256(source_bytes).hexdigest() != P242_SOURCE_SHA256:
                raise ValueError("P242 C source SHA-256 mismatch")
            source.write_bytes(source_bytes)
            radius = _export_centers(even, centers)
            subprocess.run(["cc", "-O3", "-Wall", str(source), "-o", str(binary), "-l:libmpfr.so.6", "-lgmp"],
                           check=True)
            subprocess.run([str(binary), str(centers), str(trace), "1024", str(MODES)], check=True)
            minimum_pivot, inverse_sq_point = _read_trace(trace.read_bytes(), MODES)
            hashes = {f.name: sha256(f.read_bytes()).hexdigest() for f in (centers, trace)}
    inverse_sq = _dyadic_fraction(inverse_sq_point)
    error = radius + sum(errors.values())
    scale = 10**55
    center_lower = Fraction(isqrt((scale**2 * inverse_sq.denominator - 1) // inverse_sq.numerator), scale)
    margin = center_lower - error
    if not margin > 0 or (margin + error)**2 * inverse_sq >= 1:
        raise ArithmeticError("directed MPFR replay did not certify a positive Weyl margin")
    return {
        "positive_pivots": MODES,
        "P242_source_sha256": P242_SOURCE_SHA256,
        "trace_inputs_and_outputs_sha256": hashes,
        "minimum_pivot_lower_exact": str(_dyadic_fraction(minimum_pivot)),
        "inverse_frobenius_squared_upper_exact": str(inverse_sq),
        "combined_Weyl_error_upper_exact": str(error),
        "strict_lower_margin_exact": str(margin),
    }


def _constant_floor_cap(quadrature, gram):
    """Test a common upper bound for every constant floor allowed at cut 504."""
    cap = Fraction(7, 25)
    t = Fraction(58933, 100)
    with ctx.workprec(1024):
        frequency = _ball(t)
        arch = acb(arb(1) / 4, frequency / 2).digamma().real - arb.pi().log()
        psi = arch - sum(a * (frequency * w).cos() for a, w in _prime_terms())
        if not psi < _ball(cap):
            raise ArithmeticError("the upper bound for admissible constant floors failed")
        matrix, errors = _assemble(quadrature, gram, "odd", cap)
        result = certify_center_weyl(matrix, sum(errors.values()))
        k = result.get("first_inconclusive_pivot_index")
        witness = None if k is None else _witness(matrix, k, errors)
    return {
        "admissible_beta_is_strictly_below_exact": str(cap),
        "rational_frequency_for_cap": str(t),
        "Psi4_upper_exact": str(_dyadic_fraction(psi.upper().man_exp())),
        "cap_is_used_as_a_form_upper_bound_only": True,
        "odd_validation": result,
        "negative_witness_at_the_cap": witness,
        "all_constant_floors_excluded_at_or_above_log7_over_2": (
            witness is not None and witness["negative_by_zero_extension_for_all_windows_beyond_threshold"]
        ),
    }


def postcheck(replay):
    archived, gram = _load_inputs(replay)
    old = _old_matrices(OLD_AGGREGATE)
    even = _mpfr_even(old[0], gram[0])
    upper_cap = _constant_floor_cap(old[1], gram[1])
    return {
        "schema": "NUMERIA.ES8P258.POSTCHECK.1",
        "archive_sha256": PINS,
        "even_MPFR": even,
        "GL60_endpoint_crosscheck": archived["G504_GL60_endpoint_crosscheck"],
        "GL60_condition": {s: body["condition"] for s, body in archived["sectors"].items()},
        "odd_raised_floor_exact_witness": archived["sectors"]["odd"]["negative_witness"],
        "constant_floor_cap_audit": upper_cap,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = postcheck(args.replay)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary = result["odd_raised_floor_exact_witness"]
    cap = result["constant_floor_cap_audit"]
    print("P258 SUMMARY", json.dumps({
        "even_MPFR_has_377_positive_pivots": result["even_MPFR"]["positive_pivots"] == 377,
        "even_MPFR_margin": float(Fraction(result["even_MPFR"]["strict_lower_margin_exact"])),
        "odd_rayleigh_upper_at_below_threshold_width": float(Fraction(summary["rayleigh_upper_below_threshold_exact"])),
        "odd_exact_dyadic_test_vector": summary["coordinates_exact_dyadic"],
        "all_admissible_constant_floors_excluded": cap["all_constant_floors_excluded_at_or_above_log7_over_2"],
        "odd_cap_upper_at_below_threshold_width": (
            float(Fraction(cap["negative_witness_at_the_cap"]["rayleigh_upper_below_threshold_exact"]))
            if cap["negative_witness_at_the_cap"] else None
        ),
        "condition_lower": {s: float(_dyadic_fraction(record["condition_number_lower_dyadic"]))
                            for s, record in result["GL60_condition"].items() if record is not None},
        "checked_entries": {s: r["endpoint_enclosures_contained"]
                            for s, r in result["GL60_endpoint_crosscheck"].items()},
    }))


if __name__ == "__main__":
    main()
