"""Locate the prescribed-margin obstruction before deducting the penalty B."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from .assembly import band_gram, deficit_endpoints
from .interval_ldl import dyadic


DEFAULT_HEAD = Path(__file__).with_name("prefix_replay_100.json")


def _number(pair):
    m, e = map(int, pair)
    return arb(m) * arb(2)**e


def _bounds(x):
    return [dyadic(x.lower()), dyadic(x.upper())]


def run(path=DEFAULT_HEAD):
    payload = path.read_bytes()
    data = json.loads(payload)
    if data["schema"] != "NUMERIA.P93V40.PREFIXES.1":
        raise ValueError("expected a prefix replay with exact dyadic witnesses")
    detail = [row["preconditioned_congruence"] for row in data["sectors"]]
    if any(item["status"] != "strictly_below_requested_margin" for item in detail):
        raise ValueError("the input must contain the two verified margin obstructions")
    count = max(item["unproved_point_pivot"] for item in detail)
    with ctx.workprec(768):
        ends, _ = deficit_endpoints()
        blocks = band_gram(ends, count)
        results = []
        for p, item in enumerate(detail):
            rayleigh = item["rayleigh_check"]
            witness = [_number(pair) for pair in rayleigh["witness_dyadic"]]
            if any(not x.is_exact() or not x.is_zero() for x in witness[count:]):
                raise ArithmeticError("the witness extends outside the rebuilt bands")
            if any(not x.is_exact() for x in witness[:count]):
                raise ArithmeticError("the witness must remain exact")
            vector = arb_mat(count,1,witness[:count])
            norm = (vector.transpose()*vector)[0,0]
            if not norm > 0:
                raise ArithmeticError("the zero vector is not a witness")
            lower, upper = map(_number, rayleigh["rayleigh_numerator_dyadic"])
            head_numerator = lower.union(upper)
            penalty = (vector.transpose() * (blocks[p] * vector))[0,0]
            unpenalized = (head_numerator + penalty) / norm
            if not unpenalized < arb(1)/10**60:
                raise ArithmeticError("A0 does not certify the expected margin obstruction")
            results.append({
                "parity": p,
                "support_modes": item["unproved_point_pivot"],
                "head_rayleigh_dyadic": _bounds(head_numerator/norm),
                "penalty_rayleigh_dyadic": _bounds(penalty/norm),
                "a0_rayleigh_dyadic": _bounds(unpenalized),
                "a0_strictly_below_requested_margin": True,
            })
        return {
            "schema": "NUMERIA.P93V40.OBSTRUCTION_ORIGIN.1",
            "proof_input_sha256": sha256(payload).hexdigest(),
            "target_exact": "1e-60",
            "sectors": results,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = run()
    if args.verify is not None:
        if result != json.loads(args.verify.read_text()):
            raise ValueError("the independent penalty subtraction audit differs")
        print(f"Verified {args.verify}")
    else:
        args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
