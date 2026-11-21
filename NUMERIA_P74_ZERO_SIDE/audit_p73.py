"""Replay P73 and separately verify the committed zero brackets."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from flint import acb, arb, ctx


P73_COMMIT = "e4c5cf326f8d7e05795ab5457a37deee6beb7a72"
EXPORTER = "NUMERIA_P72_ZERO_SIDE/export_zeros_flint.py"
VERIFIER = "NUMERIA_P73_ADAPTIVE/certify_adaptive.py"
ZERO_FILE = "NUMERIA_P72_ZERO_SIDE/arb_zeros_20.json"
CERTIFICATE = "NUMERIA_P73_ADAPTIVE/certificate_20_t46.6.json"
PINNED_FILES = {
    EXPORTER: "0f7ac15589d3b8032df04dfeb9fc8213a18715c1fa9d820415b94a60c4ab5e6d",
    "NUMERIA_P72_ZERO_SIDE/verify_spectral_rational.py":
        "6200a11f1cd8d97838dd90064406027fa4986d852b640aea428359523c84ce57",
    VERIFIER: "c3ce2de4d060d290f0de3de25876d3fdfbfe083b7f9eabf0b3ce827923e6bcd5",
    ZERO_FILE: "c5993b75a6db7803251655955859cd935f4381b2c0bc20fd77b9f0258bdeddab",
    CERTIFICATE: "9cca4b572b00ecb84545a5c4a34d861d96fc4d1fc7a930b9225d5727eb0d357d",
}


def _dyadic_fraction(pair):
    mantissa, exponent = map(int, pair)
    return mantissa * Fraction(2)**exponent


def _brackets(payload):
    rows = payload["zeros"]
    if len(rows) != 20 or [row["index"] for row in rows] != list(range(1, 21)):
        raise ArithmeticError("P73 replay must contain exactly the first 20 indices")
    return [
        (_dyadic_fraction(row["lo_dyad"]), _dyadic_fraction(row["hi_dyad"]))
        for row in rows
    ]


def _hardy_z(t):
    half_t = t / 2
    phase = acb(arb(1) / 4, half_t).lgamma().imag - half_t * arb.pi().log()
    return (acb(0, phase).exp() * acb(arb(1) / 2, t).zeta()).real


def _audit_hardy(brackets):
    with ctx.workprec(320):
        margin = None
        previous = arb(0)
        for index, (left, right) in enumerate(brackets, 1):
            low = arb(left.numerator) / left.denominator
            high = arb(right.numerator) / right.denominator
            if not (low.is_exact() and high.is_exact() and previous < low < high):
                raise ArithmeticError(f"invalid exact bracket at zero {index}")
            z_left, z_right = _hardy_z(low), _hardy_z(high)
            if not (z_left > 0 and z_right < 0 or z_left < 0 and z_right > 0):
                raise ArithmeticError(f"Hardy Z signs are undecided at zero {index}")
            bound = min(abs(z_left).lower(), abs(z_right).lower())
            margin = bound if margin is None else min(margin, bound)
            previous = high
        height = arb(78)
        if not previous < height or height.zeta_nzeros() != 20:
            raise ArithmeticError("total zero count does not authenticate the brackets")
        mantissa, exponent = margin.man_exp()
    return {
        "precision_bits": 320,
        "opposite_sign_brackets": 20,
        "count_height": 78,
        "count_at_height": 20,
        "minimum_absolute_hardy_endpoint_lower_dyadic": [str(mantissa), int(exponent)],
    }


def audit():
    repo = Path(__file__).resolve().parents[1]
    files = {name: (repo / name).read_bytes() for name in PINNED_FILES}
    for name, expected in PINNED_FILES.items():
        if sha256(files[name]).hexdigest() != expected:
            raise ValueError(f"P73 baseline checksum mismatch: {name}")
    original = _brackets(json.loads(files[ZERO_FILE]))
    hardy = _audit_hardy(original)

    python = [sys.executable]
    if sys.flags.optimize:
        python.append("-O")
    results = []
    with TemporaryDirectory() as scratch:
        zeros = Path(scratch) / "zeros.json"
        certificate = Path(scratch) / "certificate.json"
        for precision in (192, 256):
            subprocess.run([*python, EXPORTER, "--n", "20", "--precision", str(precision),
                            "--out", str(zeros)], cwd=repo, check=True, stdout=subprocess.DEVNULL)
            zeros_bytes = zeros.read_bytes()
            narrower = _brackets(json.loads(zeros_bytes))
            if precision == 192 and zeros_bytes != files[ZERO_FILE]:
                raise ArithmeticError("192-bit zero export differs from the pinned P72 input")
            if precision == 256 and not all(
                old_left < left < right < old_right
                for (old_left, old_right), (left, right) in zip(original, narrower, strict=True)
            ):
                raise ArithmeticError("256-bit brackets did not nest strictly inside the P72 input")

            zeros_hash = sha256(zeros_bytes).hexdigest()
            subprocess.run([*python, VERIFIER, "--zeros", str(zeros),
                            "--expected-zeros-sha256", zeros_hash,
                            "--to", "46.6", "--out", str(certificate)],
                           cwd=repo, check=True, stdout=subprocess.DEVNULL)
            certificate_bytes = certificate.read_bytes()
            if precision == 192 and certificate_bytes != files[CERTIFICATE]:
                raise ArithmeticError("192-bit rational replay differs from the pinned P73 proof")
            result = json.loads(certificate_bytes)
            if (result["cell_count"], result["zero_side_evaluations"], result["max_depth"]) != (565, 566, 12):
                raise ArithmeticError("P73 cover differs from 565 cells, 566 nodes and depth 12")
            if not Fraction(*map(int, result["min_lower_rational"])) > 0:
                raise ArithmeticError("P73 replay lower bound is nonpositive")
            results.append({
                "precision_bits": precision,
                "zeros_sha256": zeros_hash,
                "certificate_sha256": sha256(certificate_bytes).hexdigest(),
                "cell_coverage_sha256": result["cell_coverage_sha256"],
                "cell_count": result["cell_count"],
                "endpoint_evaluations": result["zero_side_evaluations"],
                "max_depth": result["max_depth"],
                "minimum_cell": result["minimum_cell"],
                "min_lower_rational": result["min_lower_rational"],
            })
    return {
        "schema": "NUMERIA.P74.P73_REPLAY.1",
        "p73_commit": P73_COMMIT,
        "start_t": "18",
        "end_t": "233/5",
        "pinned_file_sha256": PINNED_FILES,
        "zero_provenance_audit": hardy,
        "replays": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(audit(), indent=2) + "\n"
    if args.out:
        args.out.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
