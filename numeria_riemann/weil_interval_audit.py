"""Interval enclosure audit for a finite NUMERIA-RIEMANN Weil matrix.

This is a deliberately small, independent audit path for the L=0.45, N=6
prototype. It uses python-flint/Arb ball arithmetic throughout the matrix
construction, pole constraints, finite-interval integral, prime overlaps,
and a conservative analytic tail bound. It certifies positivity only if the
interval LDL pivots are all strictly positive. Otherwise it reports
"inconclusive"; it never infers a sign from mesh stability.

The archimedean integral is enclosed by interval range sums: on each cell I,
the integral of a continuous f lies in width(I)*f(I). The digamma real part
is evaluated using a positive series with an integral-test tail enclosure.
The remaining integral tail uses an elementary digamma-series bound
(see TAIL_PROOF.md), avoiding an unverified Binet remainder constant.

This certifies only the finite-dimensional constrained matrix for the stated
basis/parameters, if the interval LDL test succeeds. It is not a proof of the
global Weil criterion or of the Riemann Hypothesis.
"""
from __future__ import annotations

import argparse
import json
from decimal import Decimal, ROUND_CEILING
from pathlib import Path

from flint import arb, ctx

try:
    from .exact_mesh import exact_cell_count, validate_decimal
except ImportError:  # standalone audit execution
    from exact_mesh import exact_cell_count, validate_decimal



def _assert_ball_radius_semantics():
    """Fail closed if installed python-flint ball radius semantics differ."""
    for radius in (arb("0.1"), arb("1e-30"), arb("123456789")):
        ball = arb(0, radius)
        if not (ball.contains(radius) and ball.contains(-radius)):
            raise ArithmeticError("Arb midpoint/radius constructor failed containment test")

def _a(value):
    """Parse a user decimal as an Arb ball, avoiding binary-float constants."""
    return arb(str(value))


def _sinc(x):
    """Enclose sin(x)/x, including intervals containing zero."""
    if not x.contains(0):
        return x.sin() / x
    if not x.abs_upper() < arb("0.1"):
        raise ValueError("interval crossing sinc's removable singularity is too wide")
    total = arb(1)
    term = arb(1)
    for n in range(1, 9):
        term = term * x * x / ((2 * n) * (2 * n + 1))
        total += term if n % 2 == 0 else -term
    # For |x| <= 0.1, the alternating-series remainder after x^16/17!
    # is bounded by 0.1^18/19! < 1e-34.
    return total + arb("0 +/- 1e-34")


def _omega(r, terms):
    """Enclose Re digamma(1/4+i*r/2)-log(pi) for r>=0.

    Reflection and duplication give
      Re psi(1/4+i*r/2) = Re psi(1/2+i*r) - log(2) - pi/(2 cosh(pi*r)).
    Also Re psi(1/2+i*r) = -gamma - 2 log(2) + sum_{n>=0}
      r^2 / ((n+1/2)*((n+1/2)^2+r^2)).
    The summand decreases in n. The tail sum is between its integral
    0.5*log(1+r^2/(M+1/2)^2) and that integral plus the first omitted term.
    """
    pi = arb.pi()
    r2 = r * r
    total = -arb.const_euler() - 3 * arb(2).log() - pi.log()
    for n in range(terms):
        y = arb(n) + arb("0.5")
        total += r2 / (y * (y * y + r2))
    y = arb(terms) + arb("0.5")
    integral_tail = arb("0.5") * (1 + r2 / (y * y)).log()
    first_omitted = r2 / (y * (y * y + r2))
    # Enclose [integral_tail, integral_tail + first_omitted].
    sum_tail = integral_tail + arb(first_omitted / 2, first_omitted / 2)
    q = (-pi * r).exp()
    sech_term = pi * q / (1 + q * q)
    return total + sum_tail - sech_term


def _fourier_amplitude(r, L, k):
    """Real amplitude in F_k(r)=i*exp(-i*k*pi/2)*amplitude."""
    pi = arb.pi()
    w = arb(k) * pi / (2 * L)
    return 2 * L.sqrt() * w / (w + r) * _sinc(L * (r - w))


def _archimedean_matrix(L, N, R, step, terms, _input_R, _input_step):
    """Rigorous entrywise enclosure of (1/pi) integral_0^R Re(F_i conj(F_j))*Omega dr."""
    count = exact_cell_count(_input_R, _input_step)
    width = R / count
    pi = arb.pi()
    G = [[arb(0) for _ in range(N)] for _ in range(N)]
    for cell in range(count):
        center = (arb(2 * cell + 1) * R) / (2 * count)
        r = center + arb(0, width / 2)
        omega = _omega(r, terms)
        amplitudes = [_fourier_amplitude(r, L, k) for k in range(1, N + 1)]
        for i in range(N):
            for j in range(i, N):
                difference = (i + 1) - (j + 1)
                phase = 0 if difference % 2 else (-1) ** (difference // 2)
                if phase == 0:
                    contribution = arb(0)
                else:
                    integrand = phase * amplitudes[i] * amplitudes[j] * omega / pi
                    contribution = width * integrand
                G[i][j] += contribution
                if i != j:
                    G[j][i] += contribution
    # Analytic tail enclosure for each entry. For r >= R >= 2*w_max:
    # |F_k(r)| <= 4*w_k/(sqrt(L)*r^2), and
    # |Omega(r)| <= 8+log(3r) for r>=1; proof in TAIL_PROOF.md.
    wmax = arb(N) * pi / (2 * L)
    if not (R.lower() >= (2 * wmax).upper()):
        raise ValueError("R must be rigorously >= 2*w_max for the tail bound")
    if not (R.lower() >= 1):
        raise ValueError("R must be rigorously >= 1 for elementary tail bound")
    # Integral_R^infty [8+log(3r)]/r^4 dr
    tail_factor = (8 + (3 * R).log() + arb(1) / 3) / (3 * R ** 3)
    for i in range(N):
        wi = arb(i + 1) * pi / (2 * L)
        for j in range(i, N):
            wj = arb(j + 1) * pi / (2 * L)
            tail = 16 * wi * wj / (pi * L) * tail_factor
            G[i][j] += arb(0, tail)
            if i != j:
                G[j][i] += arb(0, tail)
    return G, count


def _von_mangoldt_prime(n):
    m = n
    p = None
    d = 2
    while d * d <= m:
        if m % d == 0:
            p = d
            while m % d == 0:
                m //= d
            break
        d += 1
    if p is None:
        return n
    return p if m == 1 else None


def _overlap(a, L, i, j):
    """Exact analytic overlap enclosure, symmetrized as in the prototype."""
    pi = arb.pi()
    wi = arb(i) * pi / (2 * L)
    wj = arb(j) * pi / (2 * L)
    b = 2 * L - a

    def directed(k, m):
        wk = arb(k) * pi / (2 * L)
        wm = arb(m) * pi / (2 * L)
        first_phase = wk * a
        if k == m:
            first = b * first_phase.cos()
        else:
            delta = wk - wm
            first = ((delta * b + first_phase).sin() - first_phase.sin()) / delta
        total = wk + wm
        second = ((total * b + first_phase).sin() - first_phase.sin()) / total
        return (first - second) / (2 * L)

    return (directed(i, j) + directed(j, i)) / 2


def _prime_matrix(L, N):
    twoL = 2 * L
    # Search ends only after Arb proves log(n) > 2L; no float exp cutoff.
    # Decimal is used only to propose a search ceiling. The Arb comparisons
    # below decide membership rigorously; avoid float overflow here.

# Limite superior para a busca de potências de primos.
   upper = (L * 2).exp().upper()
   nstop = int(float(upper)) + 2


    P = [[arb(0) for _ in range(N)] for _ in range(N)]
    active = []
    n = 2
    while True:
        logn = arb(n).log()
        if logn.lower() > twoL.upper():
            break
        p = _von_mangoldt_prime(n)
        if p is None:
            n += 1
            continue
        if logn.upper() < twoL.lower():
            active.append(n)
            coefficient = -2 * arb(p).log() / arb(n).sqrt()
            a = logn
            for i in range(1, N + 1):
                for j in range(i, N + 1):
                    value = coefficient * _overlap(a, L, i, j)
                    P[i - 1][j - 1] += value
                    if i != j:
                        P[j - 1][i - 1] += value
        elif logn.lower() <= twoL.upper():
            raise ValueError(f"prime-power cutoff is not separated from n={n}")
        n += 1
    return P, active


def _constraint_nullspace(L, N):
    """Interval-enclosed basis B for the exact two pole constraints C*B=0."""
    pi = arb.pi()
    rootL = L.sqrt()
    plus, minus = [], []
    for k in range(1, N + 1):
        w = arb(k) * pi / (2 * L)
        parity = (-1) ** k
        denom = rootL * (w * w + arb("0.25"))
        plus.append(w * ((-L / 2).exp() - parity * (L / 2).exp()) / denom)
        minus.append(w * ((L / 2).exp() - parity * (-L / 2).exp()) / denom)
    c00, c01 = plus[0], plus[1]
    c10, c11 = minus[0], minus[1]
    determinant = c00 * c11 - c01 * c10
    if determinant.contains(0):
        raise ArithmeticError("2x2 pole-constraint pivot is not separated from zero")
    dimension = N - 2
    B = [[arb(0) for _ in range(dimension)] for _ in range(N)]
    for col, j in enumerate(range(2, N)):
        B[0][col] = (-plus[j] * c11 + c01 * minus[j]) / determinant
        B[1][col] = (-c00 * minus[j] + plus[j] * c10) / determinant
        B[j][col] = arb(1)
    return B


def _project(B, W):
    N, dim = len(B), len(B[0])
    A = [[arb(0) for _ in range(dim)] for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            value = arb(0)
            for k in range(N):
                for ell in range(N):
                    value += B[k][i] * W[k][ell] * B[ell][j]
            A[i][j] = value
            A[j][i] = value
    return A


def _interval_ldl_pivots(A):
    """Interval LDL^T; positive lower bounds certify positive definiteness."""
    n = len(A)
    lower = [[arb(0) for _ in range(n)] for _ in range(n)]
    diagonal = []
    for i in range(n):
        lower[i][i] = arb(1)
        pivot = A[i][i]
        for k in range(i):
            pivot -= lower[i][k] * lower[i][k] * diagonal[k]
        diagonal.append(pivot)
        if not (pivot.lower() > 0):
            return diagonal, False
        for j in range(i + 1, n):
            value = A[j][i]
            for k in range(i):
                value -= lower[j][k] * lower[i][k] * diagonal[k]
            lower[j][i] = value / pivot
    return diagonal, True


def run(L_value=0.45, N=6, R_value=160.0, step_value=0.02, terms=128, precision=192):
    if N < 4 or terms < 16 or precision < 64:
        raise ValueError("Require N>=4, terms>=16, precision>=64 bits")
    old_precision = ctx.prec
    ctx.prec = precision
    try:
        L_input, R_input, step_input = (validate_decimal(v) for v in (L_value, R_value, step_value))
        L, R, step = _a(L_input), _a(R_input), _a(step_input)
        if not (L.lower() > 0 and R.lower() > 0 and step.lower() > 0):
            raise ValueError("L, R and step must be strictly positive")
        if not ((L * step / 2).upper() < arb("0.1")):
            raise ValueError("step too wide for sinc enclosure")
        wmax = arb(N) * arb.pi() / (2 * L)
        if not (R.lower() >= (2 * wmax).upper()):
            raise ValueError("R must be >= 2*w_max to use the analytic tail bound")
        _assert_ball_radius_semantics()
        G, cells = _archimedean_matrix(L, N, R, step, terms, R_input, step_input)
        P, active = _prime_matrix(L, N)
        W = [[G[i][j] + P[i][j] for j in range(N)] for i in range(N)]
        B = _constraint_nullspace(L, N)
        projected = _project(B, W)
        pivots, certified = _interval_ldl_pivots(projected)
        report = {
            "project": "NUMERIA-RIEMANN",
            "method": "Arb interval range-sum quadrature + analytic tail + interval LDL",
            "status": "certified_positive_definite_finite_matrix" if certified else "inconclusive_interval_ldl",
            "claim_scope": "finite constrained sine-basis matrix only; not global Weil positivity or RH",
            "parameters": {"L": str(L), "N": N, "R": str(R), "requested_step": str(step),
                           "quadrature_cells": cells, "digamma_series_terms": terms,
                           "arb_precision_bits": precision},
            "active_prime_powers": active,
            "projected_dimension": N - 2,
            "interval_ldl_pivots": [str(p) for p in pivots],
            "interval_ldl_pivot_lower_bounds": [p.lower().str(40, radius=False) for p in pivots],
            "archimedean_entry_enclosures": [[str(G[i][j]) for j in range(N)] for i in range(N)],
            "projected_matrix_enclosures": [[str(projected[i][j]) for j in range(N - 2)] for i in range(N - 2)],
            "limitations": [
                "The finite integral enclosure uses Arb interval range sums, not mesh-difference heuristics.",
                "The infinite archimedean tail is enclosed using the documented analytic bound.",
                "The pole constraints and prime-overlap integrals are enclosed with Arb arithmetic.",
                "A positive interval LDL result certifies only this finite-dimensional constrained matrix.",
                "A failed LDL test is inconclusive and is not evidence of a negative eigenvalue.",
                "No result here proves the global Weil criterion or the Riemann Hypothesis."
            ]
        }
        return report
    finally:
        ctx.prec = old_precision


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=str, default="0.45")
    parser.add_argument("--N", type=int, default=6)
    parser.add_argument("--R", type=str, default="160")
    parser.add_argument("--step", type=str, default="0.02")
    parser.add_argument("--series-terms", type=int, default=128)
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--out", type=Path, default=Path("weil-interval-audit.json"))
    parser.add_argument("--require-positive", action="store_true",
                        help="exit nonzero unless interval LDL certifies positive definiteness")
    args = parser.parse_args()
    report = run(args.L, args.N, args.R, args.step, args.series_terms, args.precision)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in (
        "status", "parameters", "active_prime_powers", "projected_dimension",
        "interval_ldl_pivot_lower_bounds", "limitations"
    )}, indent=2))
    if args.require_positive and report["status"] != "certified_positive_definite_finite_matrix":
        raise SystemExit("interval LDL did not certify the reference finite matrix")


if __name__ == "__main__":
    main() 
