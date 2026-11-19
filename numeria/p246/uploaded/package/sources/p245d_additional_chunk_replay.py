#!/usr/bin/env python3
"""P245D: independent fixed-point replay of selected entries on P237 chunks.

Recomputes all 60-node GL60 contributions for 13 matrix pairs, both parities,
and one complete 25-panel chunk. It compares the resulting enclosures against
pinned 900-bit P237 Arb chunk records. It does NOT recompute the full triangle
or the analytic quadrature remainder.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
P237_DIR = HERE.parents[1] / "p237"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(P237_DIR) not in sys.path:
    sys.path.insert(0, str(P237_DIR))

from p245b_certified_gl60 import (
    I, iv, S, PI, LN2, PREC, log_interval, certified_nodes, re_digamma,
    sin_cos, spherical_j_small, double_factorial_odd,
    source_ball_to_fraction,
)
from p237_core import N, UPPER, validate_headers

LQ = Q(161, 200)
PAIRS = ((0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (2, 2),
         (3, 3), (0, 16), (16, 16), (0, 376), (1, 376), (376, 376))
INDICES = sorted({index for pair in PAIRS for index in pair})
SUPPORTED_STARTS = (25, 500, 975)
PINNED_HASHES = {
    (25, "even"): "3903e507a9b7c492cdab22b8a290481d43593d69d9ca0f8fde9fb80cb6d73da8",
    (25, "odd"): "ff0f5256cd92893ee2b27c2753b9c97e8a17c43e55c92c8c76c908beb899eeda",
    (500, "even"): "b448bcaad18577b10a056970c0e277873f316b9b8ff35dffb6f09e4547505d71",
    (500, "odd"): "0380e8986f374d00fcf97f6b4d851cffbd596db89a518ae71b7eb983bf484607",
    (975, "even"): "7a7d4f3290c0ff4a14a728af37af322bc723cb73bd9c92f5f012b12d7c49ed2a",
    (975, "odd"): "d26488fdf2af0b26fa3bd7c0e7b767824c202432265b13e697979143cdc344be",
}


def triangle_index(i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return i * N - i * (i - 1) // 2 + (j - i)


def bessel_series(x: I, n: int, x_bound: Q = Q(21), terms: int = 140) -> I:
    """Entire power series for j_n(x), with exact rational tail majorant."""
    x = iv(x)
    if x.l <= 0 or x.u > x_bound * S:
        raise ValueError("series Bessel domain violated")
    x2 = x * x
    term = iv(1)
    total = term
    for k in range(1, terms):
        term = -term * x2 / Q(2 * k * (2 * n + 2 * k + 1))
        total = total + term
    q = Q(1)
    for k in range(1, terms + 1):
        q *= x_bound * x_bound / Q(2 * k * (2 * n + 2 * k + 1))
    ratio = x_bound * x_bound / Q(2 * (terms + 1) * (2 * n + 2 * (terms + 1) + 1))
    if ratio >= 1:
        raise ArithmeticError("Bessel series tail ratio is not below one")
    tail = q / (1 - ratio)
    if tail >= Q(1, S):
        raise ArithmeticError("Bessel series tail exceeds one fixed-point unit")
    total = total + I.endpoints(-tail, tail)
    return (x ** n) * total / double_factorial_odd(n)


def forward_all(x: I, max_degree: int = 33) -> list[I]:
    """Stable forward spherical-Bessel recurrence for x>30 and n<=33."""
    if x.l <= 30 * S:
        raise ValueError("forward recurrence requires x>30")
    small = spherical_j_small(x)
    values = [small[0], small[1]]
    for n in range(1, max_degree):
        values.append((2 * n + 1) * values[n] / x - values[n - 1])
    return values


def high_mode_bound(panel: int, parity: int) -> Q:
    n = 2 * 376 + parity
    x_max = Q(161 * (panel + 1), 400)
    return Q(50) * x_max**n / double_factorial_odd(n)


def load_reference(reference_dir: Path, sector: str, start: int, stop: int):
    filename = f"P237_CHUNK_{sector}_{start:04d}_{stop:04d}.json.gz"
    path = reference_dir / filename
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    expected = PINNED_HASHES[(start, sector)]
    if digest != expected:
        raise ArithmeticError(f"reference SHA-256 mismatch: {filename}")
    obj = json.loads(gzip.decompress(data))
    entries = validate_headers(obj, sector, start, stop, baseline=False)
    return obj, entries, digest, filename


def evaluate_chunk(reference_dir: Path, roots_file: Path, start: int):
    if start not in SUPPORTED_STARTS:
        raise ValueError(f"supported starts are {SUPPORTED_STARTS}")
    stop = start + 25
    refs = {}
    source_meta = {}
    for parity, sector in ((0, "even"), (1, "odd")):
        obj, entries, digest, filename = load_reference(reference_dir, sector, start, stop)
        refs[parity] = entries
        source_meta[sector] = {"filename": filename, "sha256": digest,
                               "panel_range": obj["panel_range"],
                               "upper_triangle_entries": len(entries)}

    roots = certified_nodes(roots_file)
    log_pi = log_interval(PI)
    logs = [LN2, log_interval(3), 2 * LN2, log_interval(5)]
    coeff = [2 * LN2 / iv(2).sqrt(),
             2 * logs[1] / iv(3).sqrt(),
             LN2,
             2 * logs[3] / iv(5).sqrt()]
    beta = log_interval(iv(504) / (2 * PI)) - iv(Q(1, 504)) - sum(coeff, iv(0))
    L = iv(LQ)
    norms = {
        parity: {i: (2 * L * (2 * (2 * i + parity) + 1)).sqrt()
                 for i in INDICES if i != 376}
        for parity in (0, 1)
    }
    sums = {parity: {pair: iv(0) for pair in PAIRS} for parity in (0, 1)}
    method = "power_series_140_terms" if start == 25 else "forward_recurrence_to_degree_33"

    for panel in range(start, stop):
        high = {parity: high_mode_bound(panel, parity) for parity in (0, 1)}
        for root, weight, _ in roots:
            t = (2 * panel + 1 + root) / 4
            x = L * t
            dig = re_digamma(t / 2)
            psi = dig - log_pi - sum((a * sin_cos(t * ln)[1]
                                      for a, ln in zip(coeff, logs)), iv(0))
            node_weight = weight * (psi - beta) / (4 * PI)

            needed_degrees = tuple(range(8)) + (32, 33)
            if start == 25:
                jvalues = {n: bessel_series(x, n) for n in needed_degrees}
            else:
                all_j = forward_all(x, 33)
                jvalues = {n: all_j[n] for n in needed_degrees}

            u = {}
            for i in INDICES:
                if i == 376:
                    u[i] = {parity: I.endpoints(-high[parity], high[parity])
                            for parity in (0, 1)}
                else:
                    for parity in (0, 1):
                        degree = 2 * i + parity
                        value = jvalues[degree]
                        if i % 2:
                            value = -value
                        u[i, parity] = norms[parity][i] * value

            for parity in (0, 1):
                for i, j in PAIRS:
                    ui = u[i] [parity] if i == 376 else u[i, parity]
                    uj = u[j] [parity] if j == 376 else u[j, parity]
                    sums[parity][i, j] = sums[parity][i, j] + node_weight * ui * uj
        if (panel - start + 1) % 5 == 0:
            print(f"P245D replay {start}..{stop}: panel {panel + 1}/{stop}", flush=True)

    comparisons = []
    for parity, sector in ((0, "even"), (1, "odd")):
        entry_list = refs[parity]
        records = []
        for pair in PAIRS:
            index = triangle_index(*pair)
            ii, jj, ball = entry_list[index]
            if [ii, jj] != list(pair):
                raise ArithmeticError(f"reference triangle order mismatch {sector} {pair}")
            lo, hi = source_ball_to_fraction(ball)
            replay = sums[parity][pair]
            contained = replay.contain_fraction(lo) and replay.contain_fraction(hi)
            if not contained:
                raise ArithmeticError(f"independent interval misses Arb source {sector} {pair}")
            records.append({
                "pair": list(pair),
                "Legendre_degrees": [2 * pair[0] + parity, 2 * pair[1] + parity],
                "replay_interval": replay.desc(),
                "source_Arb_ball_dyadic": ball,
                "source_interval_lower_fraction": str(lo),
                "source_interval_upper_fraction": str(hi),
                "source_Arb_interval_contained": True,
            })
        comparisons.append({"sector": sector, "source": source_meta[sector],
                            "selected_entry_comparisons": records})

    return {
        "schema": "NUMERIA.ES8P245D.SELECTED_ADDITIONAL_P237_CHUNK_REPLAY.1",
        "parameters": {"L": "161/200", "Tsharp": 504,
                       "panel_range": [start, stop], "panels": 25,
                       "GL60_nodes_per_panel": 60, "fixed_point_bits": PREC,
                       "modes_per_parity": 377, "selected_pairs_per_sector": len(PAIRS)},
        "root_brackets_sha256": hashlib.sha256(roots_file.read_bytes()).hexdigest(),
        "method": method,
        "reference_chunks": source_meta,
        "panel_pair_evaluations": 25 * 60 * len(PAIRS) * 2,
        "verified_chunk_sum_enclosures": 2 * len(PAIRS),
        "comparisons": comparisons,
        "full_377_by_377_triangle_recomputed": False,
        "analytic_GL60_remainder_added": False,
        "Weil_Zhu_identity_verified": False,
        "RH_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--references", type=Path, required=True,
                        help="directory with the pinned P237 900-bit chunk JSON.gz files")
    parser.add_argument("--roots", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True, choices=SUPPORTED_STARTS)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = evaluate_chunk(args.references, args.roots, args.start)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print("P245D: independently enclosed", report["verified_chunk_sum_enclosures"],
          "selected chunk sums across", report["parameters"]["panel_range"], flush=True)


if __name__ == "__main__":
    main()
