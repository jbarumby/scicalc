"""Arb Fourier–Legendre head with directed quadrature and exact band primitives."""

import csv
from fractions import Fraction
import json
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

from .bessel_energy import bessels_at_point


ARCHIVE = Path(__file__).resolve().parents[1] / "NUMERIA_P93V39"
WIDTH = Fraction(
    "69314718055994530941723212145817656807550013436025525412068000949339362197/"
    "50000000000000000000000000000000000000000000000000000000000000000000000000"
)
T = 4500
LOCALIZATION = 600
PRIMES = ((2, 2), (3, 3), (4, 2), (5, 5), (7, 7), (8, 2), (9, 3), (11, 11), (13, 13))
GL_ORDER = 100


def _ball(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def prime16_majorant(modes, parity):
    """Dominate the border term on a finite, parity-preserving Legendre span."""
    a = _ball(WIDTH)
    delta = 2 * a - arb(16).log()
    if not (delta > 0 and 2 * a < arb(17).log()):
        raise ArithmeticError("the width must cross only the threshold at 16")
    kappa = arb(2).log() / 4
    orders = range(parity, 2 * modes, 2)
    endpoints = [((2 * n + 1) / (2 * a)).sqrt() for n in orders]
    derivative_square = sum((n**4 * (2 * n + 1) for n in orders), 0) / (2 * a**3)
    rank_one = 4 * kappa * delta
    diagonal = 4 * kappa * delta**3 * derivative_square
    result = arb_mat(modes, modes)
    for i, x in enumerate(endpoints):
        for j in range(i + 1):
            entry = rank_one * x * endpoints[j] + (diagonal if i == j else 0)
            result[i, j] = result[j, i] = entry
    return result


def _bessels_at_point(x, degree):
    """Enclose j_0,...,j_degree with the recurrence centered at an exact point."""
    inv = 1 / x
    sine, cosine = x.sin_cos()
    if x > degree:
        values = [sine * inv, (sine * inv - cosine) * inv]
        for n in range(1, degree):
            values.append((2 * n + 1) * inv * values[-1] - values[-2])
        return values[:degree + 1]
    scale = (arb.pi() * inv / 2).sqrt()
    order = arb(2 * degree + 1) / 2
    values = [arb(0)] * (degree + 1)
    values[-1] = scale * x.bessel_j(order)
    following = scale * x.bessel_j(order + 1)
    for n in range(degree, 0, -1):
        values[n - 1], following = (2 * n + 1) * inv * values[n] - following, values[n]
    if not values[0].overlaps(sine * inv):
        raise ArithmeticError("spherical Bessel downward recurrence misses j_0")
    return values


def _bessels(x, degree):
    """Use a guarded point recurrence and then transfer the input interval."""
    midpoint = x.mid()
    input_error = arb(0, (x.rad() / 2).upper())
    if degree >= 100:
        values = bessels_at_point(midpoint, degree)
        return [value + input_error for value in values]
    bits = ctx.prec
    tolerance = arb(2)**(48-bits)
    for guard in (32, 128, 384, 896, 1920, 3968, 8064, 16256):
        with ctx.workprec(bits + guard):
            values = _bessels_at_point(midpoint, degree)
            if all(value.is_finite() and value.rad() < tolerance for value in values):
                break
    else:
        raise ArithmeticError(f"unable to enclose spherical Bessel recurrence at x={midpoint}")
    # j_n(x)=(2 i**n)^(-1) integral_{-1}^1 P_n(u) exp(i*x*u) du.
    # On the real axis |j_n'(x)| <= integral |u| du / 2 = 1/2.
    return [value + input_error for value in values]


def spectral_panel_sum(modes, *, batch_nodes=256, order=GL_ORDER):
    """Enclose the signed integral sum and a per-entry analytic GL error."""
    a, pi = _ball(WIDTH), arb.pi()
    if not a < arb(3) / 2:
        raise ArithmeticError("the analytic ellipse bounds require a < 3/2")
    terms = [(2 * arb(p).log() / arb(q).sqrt(), arb(q).log()) for q, p in PRIMES]
    if not sum((amp for amp, _ in terms), arb(0)) < 10:
        raise ArithmeticError("the analytic ellipse bounds require total prime mass < 10")
    nodes = [arb.legendre_p_root(order, j, weight=True) for j in range(order)]
    if not (all(nodes[i][0] > nodes[i + 1][0] for i in range(order - 1))
            and all(weight > 0 for _, weight in nodes)):
        raise ArithmeticError("invalid GL roots or positive weights")
    totals = [arb_mat(modes, modes), arb_mat(modes, modes)]
    values, weighted = [[], []], [[], []]
    norms = [(-1)**(n // 2) * (2 * a * (2 * n + 1)).sqrt() for n in range(2 * modes)]

    def flush():
        count = len(values[0]) // modes
        if count:
            for sector in (0, 1):
                x = arb_mat(count, modes, values[sector])
                wx = arb_mat(count, modes, weighted[sector])
                totals[sector] += x.transpose() * wx
                values[sector].clear()
                weighted[sector].clear()

    panels = [(Fraction(i, 4), Fraction(i + 1, 4)) for i in range(64)]
    panels += [(Fraction(i), Fraction(min(i + 16, T))) for i in range(16, T, 16)]
    log_pi = pi.log()
    quarter = arb(1) / 4
    for left, right in panels:
        center, radius = _ball((left + right) / 2), _ball((right - left) / 2)
        for node, weight in nodes:
            t = center + radius * node
            j = _bessels(a * t, 2 * modes - 1)
            symbol = acb(quarter, t / 2).digamma().real - log_pi - 2
            symbol -= sum((amp * (t * shift).cos() for amp, shift in terms), arb(0))
            scalar = radius * weight * symbol / pi
            for sector in (0, 1):
                row = [norms[n] * j[n] for n in range(sector, 2 * modes, 2)]
                values[sector].extend(row)
                weighted[sector].extend(scalar * v for v in row)
            if len(values[0]) >= modes * batch_nodes:
                flush()
    flush()

    # Uniform entry envelopes on E_4 are 10^4 on the near panels and 10^41
    # on the far panels (see the proof). GL is exact through degree 2*order-1.
    remainder = Fraction(16 * (16 * 10**4 + (T - 16) * 10**41), 3 * 4**(2 * order))
    return totals, remainder, len(panels)


def _primitive_data(x, count):
    j = _bessels(x, count)
    derivatives = [n * j[n] / x - j[n + 1] for n in range(count)]
    diagonal = [(2 * x).si() - x.sin()**2 / x]
    for n in range(1, count):
        diagonal.append(((2 * n - 1) * diagonal[-1] - x * (j[n - 1]**2 + j[n]**2)) / (2 * n + 1))
    return j[:count], derivatives, diagonal


def deficit_endpoints(archive=ARCHIVE):
    """Collect signed endpoint weights after checking the component rule."""
    archive = Path(archive)
    with (archive / "deficit_2000000.csv").open(newline="", encoding="utf-8") as file:
        reader = csv.reader(file)
        if next(reader) != ["left16", "right16", "deficit_micro"]:
            raise ValueError("unsupported deficit table header")
        bands = [tuple(map(int, row)) for row in reader]
    components = []
    previous_right = 16 * T
    for left, right, deficit in bands:
        if right <= 16 * T:
            continue
        if not previous_right <= left < right or not deficit > 0:
            raise ValueError("invalid ordered deficit interval above the cutoff")
        if not components or left != previous_right:
            components.append([])
        components[-1].append((left, right, deficit))
        previous_right = right
    factors = json.loads((archive / "hybrid39.json").read_text(encoding="utf-8"))["ks"]
    totals = {}
    for component, factor in zip(components, factors, strict=True):
        expected = 14 if sum((right - left) * d for left, right, d in component) >= 16_000_000 else 2
        if type(factor) is not int or factor != expected:
            raise ValueError("the multiplier disagrees with the exact component mass")
        for left, right, d in component:
            totals[left] = totals.get(left, 0) - factor * d
            totals[right] = totals.get(right, 0) + factor * d
    return [(Fraction(t, 16), Fraction(w, 1_000_000)) for t, w in sorted(totals.items()) if w], {
        "component_count": len(components), "bands": sum(map(len, components)),
    }


def band_gram(endpoints, modes, *, batch_nodes=256):
    """Integrate each pair over signed bands via the spherical-Bessel Wronskian."""
    a = _ball(WIDTH)
    norm = [(-1)**(n // 2) * arb(2 * n + 1).sqrt() for n in range(2 * modes)]
    cross = [arb_mat(modes, modes), arb_mat(modes, modes)]
    diagonal = [arb(0) for _ in norm]
    left, right = [[], []], [[], []]
    pi_scale = 2 / arb.pi()

    def flush():
        rows = len(left[0]) // modes
        if rows:
            for p in (0, 1):
                x = arb_mat(rows, modes, left[p])
                y = arb_mat(rows, modes, right[p])
                cross[p] += x.transpose() * y
                left[p].clear()
                right[p].clear()

    for time, weight in endpoints:
        x = a * _ball(time)
        j, derivative, d = _primitive_data(x, 2 * modes)
        w = _ball(weight)
        for n in range(2 * modes):
            diagonal[n] += w * d[n]
            left[n % 2].append(norm[n] * j[n])
            right[n % 2].append(w * x**2 * norm[n] * derivative[n])
        if len(left[0]) >= modes * batch_nodes:
            flush()
    flush()

    matrices = [arb_mat(modes, modes), arb_mat(modes, modes)]
    for p in (0, 1):
        for i in range(modes):
            n = 2 * i + p
            matrices[p][i, i] = pi_scale * (2 * n + 1) * diagonal[n]
            for j in range(i):
                m = 2 * j + p
                value = pi_scale * (cross[p][j, i] - cross[p][i, j]) / (n * (n + 1) - m * (m + 1))
                matrices[p][i, j] = matrices[p][j, i] = value
    return matrices


def assemble(modes, bits=384, archive=ARCHIVE):
    """Enclose full A0 and G_T with B supported exactly on the first 600 modes."""
    if type(modes) is not int or modes < 1 or type(bits) is not int or bits < 256:
        raise ValueError("require positive modes and at least 256-bit precision")
    with ctx.workprec(bits):
        order = GL_ORDER if bits < 512 else 140
        weighted, remainder, panels = spectral_panel_sum(modes, order=order)
        endpoints, metadata = deficit_endpoints(archive)
        b_blocks = band_gram(endpoints, min(modes, LOCALIZATION))
        g_blocks = band_gram([(T, 1)], modes)
        a = _ball(WIDTH)
        pole_argument = a / 2
        pole_scale = (arb.pi() / (2 * pole_argument)).sqrt()
        error = arb(0, _ball(remainder).upper())
        reports = []
        for p in (0, 1):
            a0, b, g = weighted[p], b_blocks[p], g_blocks[p]
            poles = [pole_scale * (2 * a * (2 * n + 1)).sqrt()
                     * pole_argument.bessel_i(arb(2 * n + 1) / 2)
                     for n in range(p, 2 * modes, 2)]
            head = arb_mat(modes, modes)
            for i in range(modes):
                for j in range(i + 1):
                    value = a0[i, j].intersection(a0[j, i])
                    value += 2 * (-1)**p * poles[i] * poles[j] + error
                    if i == j:
                        value += 2
                    a0[i, j] = a0[j, i] = value
                    # Without harmonic correction, use this finite block for
                    # N<=600. For later modes the caller supplies Q penalty.
                    if i < b.nrows():
                        value -= b[i, j]
                    if not value.is_finite():
                        raise ArithmeticError(f"nonfinite enclosure at parity {p}, entry ({i}, {j})")
                    head[i, j] = head[j, i] = value
            reports.append({
                "a0": a0, "b": b, "g": g, "head": head,
                "remaining_penalty_after_600": modes > LOCALIZATION,
                "quadrature_entry_error_exact": str(remainder), "quadrature_panels": panels,
                "deficit_weighted_endpoints": len(endpoints), **metadata,
            })
        return tuple(reports)
