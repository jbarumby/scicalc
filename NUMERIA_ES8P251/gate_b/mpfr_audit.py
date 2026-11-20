"""Replay the Gate B center LDL and inverse bound in the P242 C/MPFR verifier."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from flint import ctx

from NUMERIA_ES8P249.legendre_tail import HALF_WIDTH, four_prime_model_half_width

from .finite_head import load_aggregate
from .finite_head_certify import assemble_finite_head
from .finite_head_quadrature import CUTOFF, MODES, SECTORS


PRECISION_BITS = 1024
P242_SOURCE_SHA256 = "3fe68cac0ba12d72b1036d42b4b3c30933c3c2db2f3f520edbf846071d217009"
P242_SOURCE = Path(__file__).resolve().parents[2] / "numeria" / "p242" / "p242_mpfr_interval_ldl.c"
SCHEMA = "NUMERIA.ES8P251.GATE_B.INDEPENDENT_MPFR_WEYL.1"
_TRACE_HEADER = "# index low_mantissa low_exp high_mantissa high_exp"
_INVERSE_PREFIX = "# inverse_frobenius_squared_upper "


def _dyadic(pair):
    mantissa, exponent = map(int, pair)
    return Fraction(mantissa) * Fraction(2) ** exponent


def _read_trace(raw, dimension):
    lines = raw.decode("ascii").splitlines()
    if len(lines) != dimension + 2 or lines[0] != _TRACE_HEADER:
        raise ValueError("MPFR trace has an incorrect header or length")
    pivot_lowers = []
    for index, line in enumerate(lines[1:-1]):
        fields = line.split()
        if len(fields) != 5 or int(fields[0]) != index:
            raise ValueError(f"invalid MPFR pivot row {index}")
        lower, upper = _dyadic(fields[1:3]), _dyadic(fields[3:5])
        if not 0 < lower <= upper:
            raise ArithmeticError(f"nonpositive or inverted MPFR pivot {index}")
        pivot_lowers.append((lower, fields[1:3]))
    if not lines[-1].startswith(_INVERSE_PREFIX):
        raise ValueError("MPFR trace is missing its inverse bound")
    inverse_bound = lines[-1][len(_INVERSE_PREFIX):].split()
    if len(inverse_bound) != 2 or _dyadic(inverse_bound) <= 0:
        raise ArithmeticError("MPFR inverse Frobenius bound must be positive")
    return min(pivot_lowers)[1], inverse_bound


def _export_centers(matrix, path):
    """Export exact centers and sum original radius bounds as fractions."""
    dimension = matrix.nrows()
    radii = [Fraction(0) for _ in range(dimension)]
    with path.open("w", encoding="ascii", newline="\n") as output:
        for column in range(dimension):
            for row in range(column + 1):
                entry, mirror = matrix[row, column], matrix[column, row]
                if (not entry.is_finite() or not mirror.is_finite()
                        or entry.mid() != mirror.mid() or entry.rad() != mirror.rad()):
                    raise ValueError(f"invalid symmetric ball at ({row}, {column})")
                mantissa, exponent = entry.mid().man_exp()
                output.write(f"{row} {column} {mantissa} {exponent}\n")
                radius = _dyadic(entry.rad().upper().man_exp())
                radii[row] += radius
                if row != column:
                    radii[column] += radius
    return max(radii)


def _weyl_margin(inverse_squared, perturbation):
    """Return a rational margin m with (m + perturbation)^2 nu strictly below one."""
    scale = 10**30
    center_units = isqrt((scale**2 * inverse_squared.denominator - 1) // inverse_squared.numerator)
    margin = Fraction((center_units - perturbation * scale) // 1, scale)
    if margin <= 0 or (margin + perturbation)**2 * inverse_squared >= 1:
        raise ArithmeticError("independent Weyl margin did not certify positivity")
    return margin


def replay(path, *, compiler="cc", half_width=HALF_WIDTH, modes=MODES, cutoff=CUTOFF):
    """Audit positivity of the recorded GL60 box with independent C arithmetic."""
    half_width = four_prime_model_half_width(half_width)
    quadratures, _, input_hash = load_aggregate(path, cutoff=cutoff, half_width=half_width, modes=modes)

    source_bytes = P242_SOURCE.read_bytes()
    source_hash = sha256(source_bytes).hexdigest()
    if source_hash != P242_SOURCE_SHA256:
        raise ValueError("P242 C source SHA-256 does not match the pinned verifier")

    sectors = {}
    with TemporaryDirectory(prefix="numeria-p251-mpfr-") as directory:
        work = Path(directory)
        source = work / P242_SOURCE.name
        source.write_bytes(source_bytes)
        binary = work / "p242_mpfr_interval_ldl"
        subprocess.run(
            [compiler, "-O3", "-Wall", str(source), "-o", str(binary), "-l:libmpfr.so.6", "-lgmp"],
            check=True,
        )

        for sector, quadrature in zip(SECTORS, quadratures, strict=True):
            centers, trace = work / f"{sector}-centers.txt", work / f"{sector}-mpfr.txt"
            with ctx.workprec(PRECISION_BITS):
                assembled, metadata = assemble_finite_head(quadrature, sector, cutoff=cutoff, half_width=half_width)
                radius_row_upper = _export_centers(assembled, centers)
            subprocess.run([str(binary), str(centers), str(trace), str(PRECISION_BITS), str(modes)], check=True)
            trace_bytes = trace.read_bytes()
            minimum_pivot, inverse_bound = _read_trace(trace_bytes, modes)
            gl60_upper = metadata["gl60_error_budget"]["row_sum_bound"]
            perturbation = radius_row_upper + gl60_upper
            margin = _weyl_margin(_dyadic(inverse_bound), perturbation)
            sectors[sector] = {
                "center_dyadics_sha256": sha256(centers.read_bytes()).hexdigest(),
                "mpfr_trace_sha256": sha256(trace_bytes).hexdigest(),
                "positive_pivot_count": modes,
                "minimum_pivot_lower_dyadic": minimum_pivot,
                "inverse_frobenius_squared_upper_dyadic": inverse_bound,
                "quadrature_and_assembly_radius_row_upper_exact": str(radius_row_upper),
                "gl60_analytic_row_upper_exact": str(gl60_upper),
                "strict_margin_lower_exact": str(margin),
                "exact_weyl_inequality_verified": True,
            }

    report = {
        "schema": SCHEMA,
        "input_GL60_aggregate_sha256": input_hash,
        "p242_c_source_sha256": source_hash,
        "cutoff": cutoff,
        "modes_per_sector": modes,
        "assembly_precision_bits": PRECISION_BITS,
        "MPFR_precision_bits": PRECISION_BITS,
        "positive_pivot_count": sum(data["positive_pivot_count"] for data in sectors.values()),
        "sectors": sectors,
        "both_finite_heads_positive": True,
        "GL60_quadrature_independently_recomputed": False,
        "assumes_classical_P237_GL60_E3_budget": True,
        "Riemann_hypothesis_proved": False,
    }
    if half_width != HALF_WIDTH:
        report["L"] = str(half_width)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="aggregate GL60 matrix file")
    parser.add_argument("--out", type=Path, required=True, help="independent replay JSON")
    parser.add_argument("--cc", default="cc", help="C compiler executable")
    parser.add_argument("--half-width", type=Fraction, default=HALF_WIDTH)
    parser.add_argument("--cutoff", type=int, default=CUTOFF)
    parser.add_argument("--modes", type=int, default=MODES)
    args = parser.parse_args()
    report = replay(
        args.matrix, compiler=args.cc, half_width=args.half_width, modes=args.modes, cutoff=args.cutoff,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
