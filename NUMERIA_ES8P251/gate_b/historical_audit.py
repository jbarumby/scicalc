"""Compare fresh GL60 blocks to authenticated P237 and P245 intervals."""

import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path

from flint import ctx

from numeria.p237.p237_core import as_ball_interval, upper_indices, validate_headers

from .finite_head_quadrature import MODES, SECTORS, quadrature_chunk


ROOT = Path(__file__).resolve().parents[2]
ARCHIVED_CHUNKS = (
    (
        (0, 25),
        "numeria/p245/gate_c/P236_CHUNK25_{}.json.gz",
        ("6cfc2d3080389751b08d62dbdcdabdbd16a8be64ea5b23993dac0f27f61383f5",
         "76a7c147bf805fca4b67f65cbfd01464cf95134b30e990d3b5f24226a4de1f3d"),
    ),
    (
        (25, 50),
        "numeria/p247/references/P237_CHUNK_{}_0025_0050.json.gz",
        ("3903e507a9b7c492cdab22b8a290481d43593d69d9ca0f8fde9fb80cb6d73da8",
         "ff0f5256cd92893ee2b27c2753b9c97e8a17c43e55c92c8c76c908beb899eeda"),
    ),
)
P245_REPLAY = "numeria/p245/gate_c/P245D_SELECTED_500_525.json"
P245_REPLAY_SHA = "597ea1038dc99f0a71074f33f1be27922f2a5e15011f5e5e447f656f604e1cd2"
P245_ROOTS = "numeria/p245/package/P245B_ROOT_BRACKETS.json"
P245_ROOTS_SHA = "bf3a05df56b98e7a947eed516e25e8342a473d3c9aa01f3d194a582580742218"


def authenticated_bytes(relative, digest):
    raw = (ROOT / relative).read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise ArithmeticError(f"SHA-256 mismatch in {relative}")
    return raw


def _dyadic(value, exponent):
    return Fraction(int(value)) * Fraction(2)**int(exponent)


def _archived_interval(ball):
    lower, upper, exponent = as_ball_interval(ball)
    return _dyadic(lower, exponent), _dyadic(upper, exponent)


def _new_interval(matrix, i, j):
    value = matrix[i, j]
    return _dyadic(*value.lower().man_exp()), _dyadic(*value.upper().man_exp())


def _contains(outer, inner):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


def _direct_edges():
    comparisons = []
    for (start, stop), name, hashes in ARCHIVED_CHUNKS:
        rebuilt = quadrature_chunk(start, stop, bits=256)
        with ctx.workprec(384):
            for sector, matrix, digest in zip(SECTORS, rebuilt, hashes, strict=True):
                relative = name.format(sector)
                raw = authenticated_bytes(relative, digest)
                old = json.loads(gzip.decompress(raw))
                entries = validate_headers(old, sector, start, stop, baseline=start == 0)
                included = 0
                for (i, j), entry in zip(upper_indices(), entries, strict=True):
                    if entry[:2] != [i, j]:
                        raise ArithmeticError(f"archived index mismatch for {relative}")
                    if j < MODES:
                        if not _contains(_new_interval(matrix, i, j), _archived_interval(entry[2])):
                            raise ArithmeticError(f"reference interval not contained for {relative} at {(i, j)}")
                        included += 1
                comparisons.append({
                    "sector": sector,
                    "panel_range": [start, stop],
                    "source": relative,
                    "source_sha256": digest,
                    "archived_arb_intervals_contained": included,
                })
    return comparisons


def _fixed_point_comparison():
    old = json.loads(authenticated_bytes(P245_REPLAY, P245_REPLAY_SHA))
    authenticated_bytes(P245_ROOTS, P245_ROOTS_SHA)
    if old["root_brackets_sha256"] != P245_ROOTS_SHA:
        raise ArithmeticError("P245 root-bracket provenance mismatch")
    start, stop = old["parameters"]["panel_range"]
    if (start, stop) != (500, 525) or old["parameters"]["fixed_point_bits"] != 256:
        raise ArithmeticError("P245 replay metadata mismatch")
    rebuilt = quadrature_chunk(start, stop, bits=256)
    counts = {}
    with ctx.workprec(384):
        for sector, matrix, ref in zip(SECTORS, rebuilt, old["comparisons"], strict=True):
            if sector != ref["sector"]:
                raise ArithmeticError("P245 parity mismatch")
            included = 0
            for entry in ref["selected_entry_comparisons"]:
                i, j = entry["pair"]
                if j >= MODES:
                    continue
                reference = _archived_interval(entry["source_Arb_ball_dyadic"])
                interval = entry["replay_interval"]
                bits = interval["scale_power_of_two"]
                fixed = _dyadic(interval["lower"], -bits), _dyadic(interval["upper"], -bits)
                if bits != 256 or not _contains(fixed, reference):
                    raise ArithmeticError(f"P245 fixed-point enclosure failed at {sector} {(i, j)}")
                fresh = _new_interval(matrix, i, j)
                if not _contains(fresh, reference):
                    raise ArithmeticError(f"fresh replay excludes P245 reference at {sector} {(i, j)}")
                included += 1
            counts[sector] = included
    return {
        "panel_range": [start, stop],
        "source": P245_REPLAY,
        "source_sha256": P245_REPLAY_SHA,
        "root_brackets_sha256": P245_ROOTS_SHA,
        "archived_arb_intervals_contained_in_both_replays": counts,
    }


def run_audit():
    return {
        "schema": "NUMERIA.ES8P251.GATE_B.HISTORICAL_GL60_AUDIT.1",
        "matrix_precision_bits": 256,
        "bessel_precision_bits": 640,
        "P236_P237_complete_retained_triangles": _direct_edges(),
        "P245_selected_entries": _fixed_point_comparison(),
        "continuous_integral_remainder_part_of_comparison": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(run_audit(), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
