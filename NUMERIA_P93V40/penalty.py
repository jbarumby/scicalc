"""Reconstruct the full leading penalty without computing the low-band head."""

import argparse
import json
from pathlib import Path

from flint import ctx

from .assembly import LOCALIZATION, WIDTH, band_gram, deficit_endpoints
from .interval_ldl import dyadic
from .prefix import _head_digest


def run(modes=LOCALIZATION, bits=384):
    if type(modes) is not int or not 0 < modes <= LOCALIZATION:
        raise ValueError("the penalty is supported in at most 600 modes")
    if type(bits) is not int or bits < 256:
        raise ValueError("use at least 256-bit working precision")
    with ctx.workprec(bits):
        endpoints, metadata = deficit_endpoints()
        blocks = band_gram(endpoints, modes)
        sectors = []
        for parity, matrix in enumerate(blocks):
            minimum = min(matrix[i,i].lower() for i in range(modes))
            radius = max(matrix[i,j].rad().upper() for i in range(modes) for j in range(i+1))
            if not minimum > 0 or not radius.is_finite():
                raise ArithmeticError("the penalty enclosure is not accurate enough")
            sectors.append({
                "parity": parity,
                "ordered_lower_triangle_sha256": _head_digest(matrix),
                "minimum_diagonal_lower_dyadic": dyadic(minimum),
                "entry_radius_upper_dyadic": dyadic(radius),
            })
        return {
            "schema": "NUMERIA.P93V40.PENALTY.1",
            "width_exact": str(WIDTH),
            "modes": modes,
            "working_precision_bits": bits,
            "signed_endpoints": len(endpoints),
            **metadata,
            "sectors": sectors,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, default=LOCALIZATION)
    parser.add_argument("--bits", type=int, default=384)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--out", type=Path)
    output.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        expected = json.loads(args.verify.read_text(encoding="utf-8"))
        if run(expected["modes"], expected["working_precision_bits"]) != expected:
            raise ValueError("reconstructed penalty differs")
        print(f"Verified {args.verify}")
    else:
        result = run(args.modes,args.bits)
        args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
