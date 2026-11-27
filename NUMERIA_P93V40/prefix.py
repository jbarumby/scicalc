"""Assemble a common maximum dimension and prove its leading prefixes in Arb."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, arb_mat, ctx

from .assembly import ARCHIVE, LOCALIZATION, WIDTH, assemble
from .harmonic import THETAS, assemble_congruence, certify_hybrid
from .interval_ldl import certify_prefixes, dyadic
from .precondition import certify_by_congruence


CHECKPOINTS = (12, 25, 50, 100, 200, 400, 600, 1000, 1500, 2000, 3309)
MARGIN = Fraction(1, 10**60)


def _head_digest(matrix):
    digest = sha256()
    for i in range(matrix.nrows()):
        row = [[dyadic(v.lower()), dyadic(v.upper())] for v in (matrix[i,j] for j in range(i+1))]
        digest.update((json.dumps(row, separators=(",", ":")) + "\n").encode("ascii"))
    return digest.hexdigest()


def run(modes=25, bits=384, margin=MARGIN):
    if type(modes) is not int or modes < 1 or type(bits) is not int or bits < 256:
        raise ValueError("require positive integral modes and at least 256 bits")
    if isinstance(margin, (float, bool)):
        raise ValueError("use an exact nonnegative margin")
    margin = Fraction(margin)
    if margin < 0:
        raise ValueError("use an exact nonnegative margin")
    constants = json.loads((ARCHIVE / "hybrid39.json").read_text(encoding="utf-8"))
    c, d, lam = (constants[key] for key in ("c", "d", "lambda"))
    sizes = sorted({modes, *(n for n in CHECKPOINTS if n <= modes)})
    records = []
    with ctx.workprec(bits):
        for parity, sector in enumerate(assemble(modes, bits)):
            head = sector["head"]
            result = {key: value for key, value in sector.items() if not isinstance(value, arb_mat)}
            if modes > LOCALIZATION:
                cap = Fraction(c)
                for i in range(LOCALIZATION, modes):
                    head[i,i] -= arb(cap.numerator) / cap.denominator
            ldl = certify_prefixes(head, sizes, margin)
            fallback = None if ldl["certified"] else certify_by_congruence(head, sizes, margin)
            covered = {item["size"] for item in ldl["prefixes"] if item["status"] == "certified"}
            if fallback:
                covered.update(item["size"] for item in fallback["prefixes"] if item["certified"])
            if modes <= LOCALIZATION:
                hybrids = [{
                    "theta_exact": str(theta),
                    "status": "certified" if modes in covered else ldl["status"],
                    "same_head_as_localized": True, "pivots_proved": ldl["pivots_proved"],
                } for theta in THETAS]
            else:
                g = sector["g"]
                hq = arb_mat([[int(i==j)-g[i,j] for j in range(LOCALIZATION, modes)]
                              for i in range(LOCALIZATION, modes)])
                data = certify_hybrid(sector["a0"], sector["b"], hq, c, d, lam, margin)
                hybrids = []
                for item in data["results"]:
                    theta = item["theta_exact"]
                    passed = item["ldl"]["certified"] or (theta == "0" and modes in covered)
                    congruence = None
                    if not passed and theta != "0":
                        weighted = assemble_congruence(sector["a0"], sector["b"], hq,
                                                       c, d, lam, theta, margin)
                        congruence = certify_by_congruence(weighted, [modes], 0)
                        passed = congruence["status"] == "certified"
                    hybrids.append({
                        "theta_exact": theta,
                        "status": "certified" if passed else item["ldl"]["status"],
                        "same_head_as_localized": theta == "0",
                        "pivots_proved": item["ldl"]["pivots_proved"],
                        "first_obstruction": item["ldl"]["first_obstruction"],
                        "preconditioned_congruence": congruence,
                    })
            records.append({
                "parity": parity,
                "head_lower_triangle_sha256": _head_digest(head),
                **result,
                **ldl,
                "verified_positive_prefixes": sorted(covered),
                "preconditioned_congruence": fallback,
                "hybrid_weights": hybrids,
            })
    return {
        "schema": "NUMERIA.P93V40.PREFIXES.1",
        "modes_assembled_per_parity": modes,
        "working_precision_bits": bits,
        "width_exact": str(WIDTH),
        "localization_modes": LOCALIZATION,
        "scope": "finite head; full crossing additionally requires the complement estimates and every needed gate",
        "sectors": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=int, default=25)
    parser.add_argument("--bits", type=int, default=384)
    parser.add_argument("--margin", type=Fraction, default=MARGIN)
    io = parser.add_mutually_exclusive_group(required=True)
    io.add_argument("--out", type=Path)
    io.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.verify is not None:
        previous = json.loads(args.verify.read_text())
        result = run(previous["modes_assembled_per_parity"], previous["working_precision_bits"],
                     previous["sectors"][0]["margin_exact"])
        if json.dumps(result, sort_keys=True) != json.dumps(previous, sort_keys=True):
            raise ValueError("prefix certificate differs from rigorous replay")
        print(f"Verified {args.verify}")
    else:
        result = run(args.modes, args.bits, args.margin)
        args.out.write_text(json.dumps(result, indent=2)+"\n")
        print([(r["parity"], r["verified_positive_prefixes"]) for r in result["sectors"]])


if __name__ == "__main__":
    main()
