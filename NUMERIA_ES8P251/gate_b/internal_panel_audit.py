"""Independent scalar cross-check for one internal P251 Gate B GL60 panel.

This is a small pilot for modes 0 and 1 in each parity sector. It compares the
production Arb-root/downward-Bessel recurrence with an independently written
scalar evaluation using authenticated rational GL60 root brackets and closed
forms for spherical Bessel degrees 0 through 3. It does not certify the full
259-mode head or the analytic quadrature remainder.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import acb, arb, ctx

from .finite_head_quadrature import quadrature_chunk
from .historical_audit import P245_ROOTS, P245_ROOTS_SHA

ROOT = Path(__file__).resolve().parents[2]
PANEL = 200
MODES = 2
SCHEMA = "NUMERIA.ES8P251.GATE_B.INTERNAL_PANEL_SCALAR_CROSSCHECK.1"


def authenticated_root_bytes(path: Path | None = None, expected_sha: str = P245_ROOTS_SHA) -> bytes:
    """Read root brackets only if their complete byte stream matches the pin."""
    path = ROOT / P245_ROOTS if path is None else Path(path)
    raw = path.read_bytes()
    actual = sha256(raw).hexdigest()
    if actual != expected_sha:
        raise ArithmeticError(f"P245 root-bracket SHA-256 mismatch: {actual}")
    return raw


def _amplitudes(z: arb, length: arb) -> tuple[arb, arb, arb, arb]:
    """Closed spherical-Bessel amplitudes for degrees 0 through 3."""
    sine, cosine = z.sin_cos()
    return (
        length.sqrt() * sine / z,
        (3 * length).sqrt() * (sine / z**2 - cosine / z),
        -(5 * length).sqrt() * ((3 / z**3 - 1 / z) * sine - 3 * cosine / z**2),
        -(7 * length).sqrt() * ((15 / z**4 - 6 / z**2) * sine + (1 / z - 15 / z**3) * cosine),
    )


def _scalar_panel(panel: int, root_document: dict, *, precision_bits: int = 384):
    """Directly sum one panel from rational root brackets and scalar formulas."""
    if type(panel) is not int or panel != PANEL:
        raise ValueError(f"this pilot is pinned to panel {PANEL}")
    length = arb(161) / 100  # length = 2L; Bessel argument is L*t.
    even = [[arb(0) for _ in range(MODES)] for _ in range(MODES)]
    odd = [[arb(0) for _ in range(MODES)] for _ in range(MODES)]
    with ctx.workprec(precision_bits):
        denominator = arb(2) ** int(root_document["root_denominator_power_of_two"])
        degree = int(root_document["root_degree"])
        roots = root_document["roots"]
        if degree != 60 or len(roots) != 60:
            raise ArithmeticError("root-bracket document is not a complete GL60 set")
        for index, root in enumerate(roots):
            if root.get("index") != index:
                raise ArithmeticError("root-bracket indices are not canonical")
            left = arb(int(root["lower_numerator"])) / denominator
            right = arb(int(root["upper_numerator"])) / denominator
            if not left < right or not (right - left) <= 1 / denominator:
                raise ArithmeticError("invalid rational root bracket")
            relative_node = (left + right) / 2 + arb(0, (right - left) / 2)
            # At a root of P_60, (1-x^2)P'_60(x)=60P_59(x).
            weight = 2 * (1 - relative_node**2) / (
                degree * relative_node.legendre_p(degree - 1)
            ) ** 2
            t = (2 * panel + 1 + relative_node) / 4
            amplitudes = _amplitudes(length * t / 2, length)
            # Algebraically rewrite Psi-beta using 1-cos(x)=2sin^2(x/2).
            multiplier = acb(arb(1) / 4, t / 2).digamma().real - arb(252).log() + arb(1) / 504
            for power, prime in ((2, 2), (3, 3), (4, 2), (5, 5)):
                phase = t * arb(power).log() / 2
                multiplier += 4 * arb(prime).log() / arb(power).sqrt() * phase.sin() ** 2
            scaled_weight = weight * multiplier / (4 * arb.pi())
            for parity, matrix in enumerate((even, odd)):
                for i in range(MODES):
                    for j in range(MODES):
                        matrix[i][j] += scaled_weight * amplitudes[2 * i + parity] * amplitudes[2 * j + parity]
    return even, odd


def _endpoint(value) -> list[str]:
    mantissa, exponent = value.man_exp()
    return [str(mantissa), str(exponent)]


def _bounds(value) -> list[list[str]]:
    return [_endpoint(value.lower()), _endpoint(value.upper())]


def _fraction(pair: list[str]) -> Fraction:
    mantissa, exponent = map(int, pair)
    return Fraction(mantissa) * Fraction(2) ** exponent


def _overlap(left, right) -> bool:
    left_lo, left_hi = map(_fraction, left)
    right_lo, right_hi = map(_fraction, right)
    return max(left_lo, right_lo) <= min(left_hi, right_hi)


def build_report(*, panel: int = PANEL, root_path: Path | None = None) -> dict:
    if type(panel) is not int or panel != PANEL:
        raise ValueError(f"this pilot is pinned to panel {PANEL}")
    root_raw = authenticated_root_bytes(root_path)
    root_document = json.loads(root_raw)
    reference = _scalar_panel(panel, root_document)
    production_even, production_odd = quadrature_chunk(
        panel, panel + 1, modes=MODES, bits=256, bessel_bits=640, threads=1, batch_panels=1,
    )
    sectors = []
    for sector, expected, obtained in zip(
        ("even", "odd"), reference, (production_even, production_odd), strict=True,
    ):
        entries = []
        for i in range(MODES):
            for j in range(i, MODES):
                direct_bounds = _bounds(expected[i][j])
                production_bounds = _bounds(obtained[i, j])
                overlaps = _overlap(direct_bounds, production_bounds)
                if not overlaps:
                    raise ArithmeticError(f"independent scalar interval misses production at {sector}[{i},{j}]")
                entries.append({
                    "pair": [i, j],
                    "Legendre_degrees": [2 * i + (sector == "odd"), 2 * j + (sector == "odd")],
                    "scalar_reference_dyadic_bounds": direct_bounds,
                    "production_dyadic_bounds": production_bounds,
                    "intervals_overlap": True,
                })
        sectors.append({"sector": sector, "entries": entries})
    report = {
        "schema": SCHEMA,
        "panel_range": [panel, panel + 1],
        "modes_per_parity": MODES,
        "GL60_nodes": 60,
        "matrix_precision_bits": 256,
        "scalar_reference_precision_bits": 384,
        "root_brackets_sha256": sha256(root_raw).hexdigest(),
        "reference_method": "authenticated rational P245 GL60 root brackets; direct scalar formulas for j0..j3; beta identity rewritten with sin^2",
        "production_method": "P251 quadrature_chunk; Arb-native GL60 roots and downward spherical-Bessel recurrence",
        "sectors": sectors,
        "all_selected_intervals_overlap": True,
        "beta_identity_poles_and_analytic_remainder_added": False,
        "full_259_mode_head_recomputed": False,
        "GL60_analytic_remainder_certified": False,
        "finite_head_positivity_recertified": False,
        "Weil_Zhu_identity_certified": False,
        "Riemann_hypothesis_proved": False,
    }
    verify_report(report)
    return report


def verify_report(report: dict, *, expected_root_sha: str = P245_ROOTS_SHA) -> None:
    """Fail closed on schema, pinned panel, entry indices, or interval mismatch."""
    if report.get("schema") != SCHEMA or report.get("panel_range") != [PANEL, PANEL + 1]:
        raise ValueError("wrong internal-panel certificate schema or panel range")
    if report.get("root_brackets_sha256") != expected_root_sha:
        raise ArithmeticError("internal-panel certificate root hash mismatch")
    if report.get("modes_per_parity") != MODES or report.get("GL60_nodes") != 60:
        raise ValueError("internal-panel certificate dimensions mismatch")
    if report.get("sectors") is None or [item.get("sector") for item in report["sectors"]] != ["even", "odd"]:
        raise ValueError("expected exactly the even and odd sectors")
    seen = set()
    for sector in report["sectors"]:
        if len(sector.get("entries", [])) != 3:
            raise ValueError("expected all three upper-triangle entries for two modes")
        for entry in sector["entries"]:
            pair = entry.get("pair")
            if pair not in ([0, 0], [0, 1], [1, 1]):
                raise ValueError("unexpected matrix entry")
            key = sector["sector"], tuple(pair)
            if key in seen:
                raise ValueError("duplicate matrix entry")
            seen.add(key)
            left, right = entry.get("scalar_reference_dyadic_bounds"), entry.get("production_dyadic_bounds")
            if (not isinstance(left, list) or len(left) != 2 or not isinstance(right, list) or len(right) != 2
                    or not _fraction(left[0]) <= _fraction(left[1])
                    or not _fraction(right[0]) <= _fraction(right[1])
                    or not _overlap(left, right) or entry.get("intervals_overlap") is not True):
                raise ArithmeticError("invalid or nonoverlapping panel intervals")
    if len(seen) != 6 or report.get("all_selected_intervals_overlap") is not True:
        raise ValueError("incomplete internal-panel certificate")
    for field in (
        "beta_identity_poles_and_analytic_remainder_added",
        "full_259_mode_head_recomputed",
        "GL60_analytic_remainder_certified",
        "finite_head_positivity_recertified",
        "Weil_Zhu_identity_certified",
        "Riemann_hypothesis_proved",
    ):
        if report.get(field) is not False:
            raise ValueError(f"scope flag {field} must remain false")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=int, default=PANEL)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(panel=args.panel)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("P251 internal-panel scalar audit: 6/6 selected entries overlap", flush=True)


if __name__ == "__main__":
    main()
