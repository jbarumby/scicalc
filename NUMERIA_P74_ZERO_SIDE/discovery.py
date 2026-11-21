"""Compare floating-point scouting strategies at fixed evaluation budgets."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import random


ZERO_FILE = Path(__file__).resolve().parents[1] / "NUMERIA_P72_ZERO_SIDE" / "arb_zeros_20.json"
ZERO_SHA256 = "c5993b75a6db7803251655955859cd935f4381b2c0bc20fd77b9f0258bdeddab"
HEIGHT = 3 * 10**12
DEFAULT_BUDGETS = (256, 512)
DEFAULT_SEEDS = (73, 74, 101)


def _zero_midpoints():
    data = ZERO_FILE.read_bytes()
    if sha256(data).hexdigest() != ZERO_SHA256:
        raise ValueError("the P72 zero file does not match its pinned SHA-256")
    result = []
    for item in json.loads(data)["zeros"]:
        lo, hi = (int(m) * Fraction(2)**e for m, e in (item["lo_dyad"], item["hi_dyad"]))
        result.append(float((lo + hi) / 2))
    return result


def _coarse_tail_mass():
    log_height = math.log(HEIGHT)
    error = 0.10076 * log_height + 0.24460 * math.log(log_height) + 8.08344
    slope = 0.10076 + 0.24460 / log_height
    return ((log_height - math.log(2 * math.pi) + 1) / (2 * math.pi * HEIGHT)
            + (2 * error + slope / 2) / HEIGHT**2)


def _adaptive(measure, left, right, budget, zero_count):
    """Refine the earliest potentially difficult cell, with regular broad sampling."""
    floor_width = (right - left) / (8 * budget)
    cells = []

    def append_cell(a, b):
        left_sum, _ = measure(a)
        right_sum, right_tail = measure(b)
        width = b - a
        proxy = min(left_sum, right_sum) - zero_count * width**2 / 4 - right_tail
        cells.append((a, b, width, proxy))

    append_cell(left, right)
    for step in range(budget - 2):
        if step % 4 == 0:
            cell = min(cells, key=lambda row: (-row[2], row[0]))
        else:
            # Cells below the scouting resolution fall behind broader candidates.
            cell = min(cells, key=lambda row: (
                row[2] <= floor_width, row[3] >= 0, row[0] if row[3] < 0 else row[3],
                -row[2], row[0],
            ))
        cells.remove(cell)
        a, b, _, _ = cell
        middle = (a + b) / 2
        if not a < middle < b:
            raise ArithmeticError("floating point scouting resolution exhausted")
        append_cell(a, middle)
        append_cell(middle, b)


def _run_strategy(method, gammas, tail_mass, left, right, budget, seed=None):
    points = {}

    def measure(t):
        if t not in points:
            spectrum = math.fsum(4 * (math.sin(gamma * t / 2) / gamma)**2 for gamma in gammas)
            points[t] = spectrum, 2 * tail_mass * math.cosh(t / 2)
        return points[t]

    measure(left)
    measure(right)
    if method == "uniform":
        for i in range(1, budget - 1):
            measure(left + (right - left) * i / (budget - 1))
    elif method == "random":
        rng = random.Random(seed)
        for _ in range(8 * budget):
            if len(points) == budget:
                break
            measure(left + (right - left) * rng.random())
    elif method == "adaptive":
        _adaptive(measure, left, right, budget, len(gammas))
    else:
        raise ValueError("unknown scouting strategy")
    if len(points) != budget:
        raise ArithmeticError("scouting strategies used different evaluation budgets")

    scores = [(point, positive - loss) for point, (positive, loss) in sorted(points.items())]
    minimum_t, minimum = min(scores, key=lambda row: (row[1], row[0]))
    first_negative = next((i for i, row in enumerate(scores) if row[1] < 0), None)
    before = scores[first_negative - 1][0] if first_negative else None
    after = scores[first_negative][0] if first_negative is not None else None
    trace = "".join(f"{point:.14g},{score:.12f}\n" for point, score in scores)
    return {
        "method": method,
        "seed": seed,
        "evaluations": len(points),
        "negative_proxy_samples": sum(value < 0 for _, value in scores),
        "first_negative_proxy_t": None if after is None else format(after, ".12g"),
        "previous_nonnegative_proxy_t": None if before is None else format(before, ".12g"),
        "first_sample_bracket_width": None if before is None else format(after - before, ".8g"),
        "minimum_proxy_t": format(minimum_t, ".12g"),
        "minimum_proxy_score": format(minimum, ".10g"),
        "rounded_point_trace_sha256": sha256(trace.encode("ascii")).hexdigest(),
    }


def compare(budgets=DEFAULT_BUDGETS, seeds=DEFAULT_SEEDS, left=46, right=48):
    """Return diagnostics only. No scouting output is an input to certification."""
    if not (math.isfinite(left) and math.isfinite(right) and 0 <= left < right <= 100):
        raise ValueError("scouting endpoints must satisfy 0 <= left < right <= 100")
    if not budgets or any(type(budget) is not int or budget < 3 for budget in budgets):
        raise ValueError("each scouting budget must be an integer of at least three")
    gammas = _zero_midpoints()
    tail_mass = _coarse_tail_mass()
    strategies = [("uniform", None), ("adaptive", None), *(("random", seed) for seed in seeds)]
    return {
        "schema": "NUMERIA.P74.DISCOVERY.1",
        "scouting_interval": [left, right],
        "zero_count": len(gammas),
        "source_zero_sha256": ZERO_SHA256,
        "proof_input": False,
        "runs": [
            _run_strategy(method, gammas, tail_mass, left, right, budget, seed)
            for budget in budgets for method, seed in strategies
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", action="append", type=int, dest="budgets")
    parser.add_argument("--seed", action="append", type=int, dest="seeds")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = compare(args.budgets or DEFAULT_BUDGETS, args.seeds or DEFAULT_SEEDS)
    content = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.write_text(content, encoding="utf-8")
    else:
        print(content, end="")


if __name__ == "__main__":
    main()
