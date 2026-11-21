"""Regenerate the P74 experiment matrix and two continuous proof comparisons."""

import argparse
import json
from pathlib import Path

from flint import arb, ctx

from NUMERIA_ES8P259.suzuki import _dyadic
from .certify import certify
from .tails import TailBounds, count_tail_upper
from .zeros import verified_ordinates


COUNTS = (20, 50, 100, 200, 500)
PRECISIONS = (192, 256)
CHECKPOINTS = ("46.6", "47.48", "50", "55")


def coarse_obstruction():
    """Bound the maximum of any retained positive-line sum at t=50."""
    with ctx.workprec(192):
        ordinates, metadata = verified_ordinates(500)
        height = metadata["count_height"]
        if height is None:
            raise ArithmeticError("the obstruction calculation needs an integer count cutoff")
        total_upper = (4 * (sum((1 / z**2 for z in ordinates), arb(0))
                            + count_tail_upper(height))).upper()
        tail = TailBounds("coarse")
        tail_lower = tail.coarse(arb(50)).lower()
        if not total_upper < tail_lower:
            raise ArithmeticError("the coarse-tail obstruction was not proved")
        return {
            "zero_count_cutoff": height,
            "upper_all_line_terms_dyadic": _dyadic(total_upper),
            "coarse_tail_lower_at_50_dyadic": _dyadic(tail_lower),
            "time": 50,
        }


def run():
    graded = [
        certify(end="60", zero_count=n, precision=prec, checkpoints=CHECKPOINTS)
        for prec in PRECISIONS for n in COUNTS
    ]
    coarse = [
        certify(end="47.48", zero_count=n, precision=192, tail_profile="coarse",
                checkpoints=("46.6",))
        for n in COUNTS
    ]
    uniform = certify(zero_count=500, precision=192, method="uniform", mesh_width="1/200",
                      checkpoints=CHECKPOINTS)
    extended = certify(end="67.5", zero_count=500, precision=256,
                       checkpoints=("60",))
    return {
        "schema": "NUMERIA.P74.EXPERIMENTS.1",
        "graded_by_precision_and_count": graded,
        "coarse_by_count": coarse,
        "uniform_control": uniform,
        "graded_extension": extended,
        "coarse_obstruction": coarse_obstruction(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--out", type=Path)
    mode.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2) + "\n"
    if args.verify is not None:
        if result != args.verify.read_text(encoding="utf-8"):
            raise ValueError("experiment certificate differs from rigorous replay")
        print(f"Verified {args.verify}")
    else:
        args.out.write_text(result, encoding="utf-8")


if __name__ == "__main__":
    main()
