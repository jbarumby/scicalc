"""Arb replay of the P236 GL60 prime component and its cell projection."""

import hashlib
import json
from fractions import Fraction
from itertools import combinations_with_replacement
from numbers import Integral
from pathlib import Path

from flint import acb, arb, ctx


def _validate_panel(panel, cells, bits):
    if not isinstance(panel, Integral) or not 0 <= panel < 1008:
        raise ValueError("panel must be an integer between 0 and 1007")
    if not isinstance(cells, Integral) or cells < 1:
        raise ValueError("cells must be a positive integer")
    if not isinstance(bits, Integral) or bits < 128:
        raise ValueError("use at least 128 bits of Arb precision")


def _panel_nodes(panel):
    """Return (t, GL60 weight / pi, signed prime cosine sum) tuples."""
    quarter = arb(1) / 4
    center = arb(2 * panel + 1) / 4
    scale = quarter / arb.pi()
    prime = [
        (2 * arb(p).log() / arb(q).sqrt(), arb(q).log())
        for q, p in ((2, 2), (3, 3), (4, 2), (5, 5))
    ]
    nodes = []
    for k in range(60):
        root, weight = arb.legendre_p_root(60, k, weight=True)
        t = center + quarter * root
        oscillation = sum((a * (t * shift).cos() for a, shift in prime), arb(0))
        nodes.append((t, scale * weight, oscillation))
    prime_mass = sum((a for a, _ in prime), arb(0))
    return nodes, prime_mass


def _cell_lags(nodes, width, cells):
    # The squared sinc comes from exactly integrating the two normalized cell
    # indicators against the exponential at the same GL60 node.
    contributions = []
    for t, weight, oscillation in nodes:
        phase_step = t * width
        amplitude = width * weight * oscillation * (phase_step / 2).sinc() ** 2
        contributions.append((phase_step, amplitude))
    return tuple(
        sum((amplitude * (step * lag).cos() for step, amplitude in contributions), arb(0))
        for lag in range(cells)
    )


def panel_cell_lags(panel, cells, *, bits=160):
    """Enclose the GL60 prime-only cell Toeplitz lags for one P236 panel.

    The prime multiplier here is S(t). It enters P236's full C with a minus.
    """
    _validate_panel(panel, cells, bits)
    with ctx.workprec(int(bits)):
        nodes, _ = _panel_nodes(int(panel))
        return _cell_lags(nodes, (arb(161) / 100) / int(cells), int(cells))


def panel_cell_lags_float(panel, cells, *, bits=160):
    """Convert the Arb lag midpoints to binary64 for exploratory FFT use."""
    return tuple(float(lag.mid()) for lag in panel_cell_lags(panel, cells, bits=bits))


def _legendre_cell_coefficients(degree, half_width, cell_width, boundaries):
    # For n >= 1 use (P_(n+1)-P_(n-1))/(2n+1) as a primitive on [-1,1].
    if degree == 0:
        primitive = boundaries
    else:
        primitive = [
            (point.legendre_p(degree + 1) - point.legendre_p(degree - 1)) / (2 * degree + 1)
            for point in boundaries
        ]
    factor = half_width * (arb(2 * degree + 1) / (2 * half_width * cell_width)).sqrt()
    return tuple(factor * (right - left) for left, right in zip(primitive, primitive[1:]))


def _direct_toeplitz_action(lags, coefficients):
    # Direct Arb audit uses only the lags and the coefficient vector.
    return tuple(
        sum((lags[abs(i - j)] * value for j, value in enumerate(coefficients)), arb(0))
        for i in range(len(coefficients))
    )


def _legendre_quadrature_amplitudes(nodes, degrees, half_width):
    # On either parity sector, the common phase is 1 or i and cancels in
    # Re(F_n conjugate(F_m)). The amplitudes retain the alternating signs.
    norms = {
        n: (-1) ** (n // 2) * (2 * half_width * (2 * n + 1)).sqrt()
        for n in degrees
    }
    amplitudes = {n: [] for n in degrees}
    pi = arb.pi()
    for t, _, _ in nodes:
        z = half_width * t
        factor = (pi / (2 * z)).sqrt()
        for n in degrees:
            amplitudes[n].append(norms[n] * factor * z.bessel_j(arb(2 * n + 1) / 2))
    return amplitudes


def _full_weights(nodes, prime_mass):
    pi = arb.pi()
    t_sharp = arb(504)
    beta = (t_sharp / (2 * pi)).log() - 1 / t_sharp - prime_mass
    return [
        weight * (acb(arb(1) / 4, t / 2).digamma().real - pi.log() - oscillation - beta)
        for t, weight, oscillation in nodes
    ]


def _quadrature_entry(weights, left, right):
    return sum((weight * a * b for weight, a, b in zip(weights, left, right)), arb(0))


def replay_panel(panel, cells, degrees, *, bits=160):
    """Enclose selected same-parity Legendre and projected-cell entries.

    The P236 full-weight entry is the quadrature of C alone. The grid
    residual and its analytical bound concern the prime-only entries.
    Low degrees suit the default precision. Native Bessel replay of
    high degrees may need the 900 bits used by P236.
    """
    _validate_panel(panel, cells, bits)
    degrees = tuple(degrees)
    if not degrees or any(not isinstance(n, Integral) or not 0 <= n <= 753 for n in degrees):
        raise ValueError("choose integer P236 degrees from 0 to 753")
    if len(set(degrees)) != len(degrees) or len({n % 2 for n in degrees}) != 1:
        raise ValueError("choose distinct degrees within one parity")
    degrees = tuple(sorted(map(int, degrees)))
    panel, cells, bits = int(panel), int(cells), int(bits)

    with ctx.workprec(bits):
        half_width = arb(161) / 200
        support_length = 2 * half_width
        width = support_length / cells
        nodes, prime_mass = _panel_nodes(panel)
        lags = _cell_lags(nodes, width, cells)
        boundaries = tuple(arb(2 * j - cells) / cells for j in range(cells + 1))
        coefficients = {
            n: _legendre_cell_coefficients(n, half_width, width, boundaries)
            for n in degrees
        }
        actions = {n: _direct_toeplitz_action(lags, coefficients[n]) for n in degrees}
        amplitudes = _legendre_quadrature_amplitudes(nodes, degrees, half_width)
        prime_weights = [weight * oscillation for _, weight, oscillation in nodes]
        full_weights = _full_weights(nodes, prime_mass)
        # For one width-h cell, the covariance of two Lipschitz functions
        # integrates to at most h^3/12 times the two Lipschitz constants.
        # The centered exponential has Lipschitz constant |t|. On the
        # Fourier factors before and after projection use sqrt(2L) by
        # Cauchy–Schwarz and contraction of the orthogonal projection.
        grid_factor = (support_length * width**2 * support_length.sqrt() / 12) * sum(
            (abs(weight) * abs(t) for (t, _, _), weight in zip(nodes, prime_weights)),
            arb(0),
        )
        derivatives = {
            n: (arb(2 * n + 1) / support_length).sqrt() * n**2 / half_width
            for n in degrees
        }

        entries = []
        for n, m in combinations_with_replacement(degrees, 2):
            legendre = _quadrature_entry(prime_weights, amplitudes[n], amplitudes[m])
            full = _quadrature_entry(full_weights, amplitudes[n], amplitudes[m])
            projected = sum((a * b for a, b in zip(coefficients[n], actions[m])), arb(0))
            entries.append({
                "degrees": (n, m),
                "legendre_prime": legendre,
                "cell_projected_prime": projected,
                "grid_residual": legendre - projected,
                "grid_bound": (grid_factor * (derivatives[n] + derivatives[m])).upper(),
                "p236_full_quadrature": full,
            })
        return {"panel_index": panel, "cells": cells, "precision_bits": bits, "entries": tuple(entries)}


def _dyadic_interval(ball):
    lower = ball.lower().man_exp()
    upper = ball.upper().man_exp()
    return [[str(lower[0]), int(lower[1])], [str(upper[0]), int(upper[1])]]


def replay_panel_report(panel, cells, degrees, *, bits=160):
    """Return a JSON-safe replay with exact dyadic interval endpoints."""
    replay = replay_panel(panel, cells, degrees, bits=bits)
    with ctx.workprec(int(bits)):
        replay["entries"] = [
            {
                name: list(value) if name == "degrees" else _dyadic_interval(value)
                for name, value in entry.items()
            }
            for entry in replay["entries"]
        ]
    return replay


def _as_fraction(dyadic):
    mantissa, exponent = map(int, dyadic)
    if exponent >= 0:
        return Fraction(mantissa << exponent)
    return Fraction(mantissa, 1 << -exponent)


def _as_dyadic(value):
    return [str(value.numerator), 1 - value.denominator.bit_length()]


def _load_p236_archive(path):
    raw = Path(path).read_bytes()
    sample = json.loads(raw)
    if sample.get("sector") not in ("even", "odd"):
        raise ValueError("expected P236 even or odd sector")
    parity = int(sample["sector"] == "odd")
    _validate_panel(sample.get("panel_index"), 1, sample.get("precision_bits"))
    expected = {
        "schema": "NUMERIA.ES8P236.ARB_377_SINGLE_GL60_PANEL.1",
        "L": "161/200",
        "Tsharp": 504,
        "root_count": 60,
        "mode_count": 377,
        "Legendre_degree_max": 752 + parity,
        "all_C_entries_computed": True,
        "GL60_integral_analytic_remainder_added": False,
        "beta_and_poles_added": False,
    }
    for key, value in expected.items():
        if sample.get(key) != value:
            raise ValueError(f"P236 archive metadata mismatch at {key}")
    archive_entries = sample.get("gauss_interval_matrix_sample")
    if not archive_entries or not isinstance(archive_entries, list):
        raise ValueError("expected nonempty P236 archive samples")
    if len(archive_entries) != sample.get("sample_entries_exported"):
        raise ValueError("unexpected P236 archive entry count")

    original = {}
    for i, j, ball in archive_entries:
        if not isinstance(i, int) or not isinstance(j, int) or not 0 <= i <= j < 377:
            raise ValueError(f"invalid P236 sample indices {(i, j)}")
        degrees = 2 * i + parity, 2 * j + parity
        if degrees in original:
            raise ValueError(f"duplicate P236 sample indices {(i, j)}")
        midpoint, radius = map(_as_fraction, ball)
        if radius < 0:
            raise ValueError(f"negative P236 radius for degrees {degrees}")
        original[degrees] = midpoint - radius, midpoint + radius
    return sample, original, hashlib.sha256(raw).hexdigest()


def _compare_p236_endpoints(replay_entries, sample, original, archive_hash, bits):
    checks = []
    max_archive_width = Fraction(0)
    max_replay_width = Fraction(0)
    with ctx.workprec(bits):
        for degrees, value in replay_entries:
            if degrees not in original:
                raise ValueError(f"P236 archive has no sample for degrees {degrees}")
            if not value.is_finite():
                raise ArithmeticError(f"nonfinite replay quadrature for degrees {degrees}")
            source_lower, source_upper = original[degrees]
            lower, upper = map(_as_fraction, _dyadic_interval(value))
            max_archive_width = max(max_archive_width, source_upper - source_lower)
            max_replay_width = max(max_replay_width, upper - lower)
            checks.append({
                "degrees": list(degrees),
                "overlap": lower <= source_upper and source_lower <= upper,
                "archive_contained_in_replay": lower <= source_lower and source_upper <= upper,
                "replay_contained_in_archive": source_lower <= lower and upper <= source_upper,
            })
    return {
        "panel_index": sample["panel_index"],
        "sector": sample["sector"],
        "archive_json_sha256": archive_hash,
        "archive_github_sha": sample["github_sha"],
        "archive_precision_bits": sample["precision_bits"],
        "replay_precision_bits": bits,
        "max_archive_width": _as_dyadic(max_archive_width),
        "max_replay_width": _as_dyadic(max_replay_width),
        "sample_checks": checks,
    }


def compare_p236_archive(replay, path):
    """Compare the selected panel replay entries to a P236 archive exactly."""
    sample, original, archive_hash = _load_p236_archive(path)
    if sample["panel_index"] != replay["panel_index"]:
        raise ValueError("P236 archive panel index differs from the replay")
    entries = [(row["degrees"], row["p236_full_quadrature"]) for row in replay["entries"]]
    return _compare_p236_endpoints(
        entries, sample, original, archive_hash, replay["precision_bits"]
    )


def replay_p236_archive(path, *, bits=900):
    """Replay every sampled full C entry from a P236 panel JSON file.

    The direct scalar quadrature uses the distinct required modes and
    skips formation of the 377-by-377 Gram matrix. The precision matters
    for high-degree Bessel intervals, whose widths are also reported.
    """
    sample, original, archive_hash = _load_p236_archive(path)
    panel = sample["panel_index"]
    _validate_panel(panel, 1, bits)
    bits = int(bits)
    degrees = sorted({n for pair in original for n in pair})
    with ctx.workprec(bits):
        nodes, prime_mass = _panel_nodes(panel)
        weights = _full_weights(nodes, prime_mass)
        amplitudes = _legendre_quadrature_amplitudes(nodes, degrees, arb(161) / 200)
        entries = [
            ((n, m), _quadrature_entry(weights, amplitudes[n], amplitudes[m]))
            for n, m in sorted(original)
        ]
    return _compare_p236_endpoints(entries, sample, original, archive_hash, bits)
